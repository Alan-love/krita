/*
 *  SPDX-FileCopyrightText: 2017 Dmitry Kazakov <dimula73@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "KisPasteActionFactories.h"

#include "kis_config.h"
#include "kis_image.h"
#include "KisViewManager.h"
#include "kis_tool_proxy.h"
#include "kis_canvas2.h"
#include "kis_canvas_controller.h"
#include "kis_group_layer.h"
#include "kis_paint_device.h"
#include "kis_paint_layer.h"
#include "kis_shape_layer.h"
#include "kis_import_catcher.h"
#include "kis_clipboard.h"
#include "kis_selection.h"
#include "commands/kis_selection_commands.h"
#include "commands/kis_image_layer_add_command.h"
#include "KisTransformToolActivationCommand.h"
#include "kis_processing_applicator.h"
#include "kis_node_manager.h"

#include <KoDocumentInfo.h>
#include <KoSvgPaste.h>
#include <KoShapeController.h>
#include <KoShapeControllerBase.h>
#include <KoShapeManager.h>
#include <KoSelection.h>
#include <KoSelectedShapesProxy.h>
#include "kis_algebra_2d.h"
#include <KoShapeMoveCommand.h>
#include <KoShapeReorderCommand.h>
#include "kis_time_span.h"
#include "kis_keyframe_channel.h"
#include "kis_raster_keyframe_channel.h"
#include "kis_painter.h"
#include <KisPart.h>
#include <KisDocument.h>
#include <KisReferenceImagesLayer.h>
#include <KoShapeBackgroundCommand.h>
#include <KoShapeStrokeCommand.h>
#include <KoShapeBackground.h>
#include <KoShapeStroke.h>
#include <KisMainWindow.h>
#include <QApplication>
#include <QClipboard>

namespace {
QPointF getFittingOffset(QList<KoShape*> shapes,
                         const QPointF &shapesOffset,
                         const QRectF &documentRect,
                         const qreal fitRatio)
{
    QPointF accumulatedFitOffset;

    Q_FOREACH (KoShape *shape, shapes) {
        const QRectF bounds = shape->boundingRect();

        const QPointF center = bounds.center() + shapesOffset;

        const qreal wMargin = (0.5 - fitRatio) * bounds.width();
        const qreal hMargin = (0.5 - fitRatio) * bounds.height();
        const QRectF allowedRect = documentRect.adjusted(-wMargin, -hMargin, wMargin, hMargin);

        const QPointF fittedCenter = KisAlgebra2D::clampPoint(center, allowedRect);

        accumulatedFitOffset += fittedCenter - center;
    }

    return accumulatedFitOffset;
}

bool tryPasteShapes(KisPasteActionFactory::Flags flags, KisViewManager *view)
{
    bool result = false;

    KoSvgPaste paste;

    if (paste.hasShapes() && view->activeNode()->isEditable()) {
        KoCanvasBase *canvas = view->canvasBase();

        QSizeF fragmentSize;
        QList<KoShape*> shapes =
            paste.fetchShapes(canvas->shapeController()->documentRectInPixels(),
                              canvas->shapeController()->pixelsPerInch(), &fragmentSize);

        if (!shapes.isEmpty()) {
            KoShapeManager *shapeManager = canvas->shapeManager();
            shapeManager->selection()->deselectAll();

            // adjust z-index of the shapes so that they would be
            // pasted on the top of the stack
            QList<KoShape*> topLevelShapes = shapeManager->topLevelShapes();
            auto it = std::max_element(topLevelShapes.constBegin(), topLevelShapes.constEnd(), KoShape::compareShapeZIndex);
            if (it != topLevelShapes.constEnd()) {
                const int zIndexOffset = (*it)->zIndex();

                std::stable_sort(shapes.begin(), shapes.end(), KoShape::compareShapeZIndex);

                QList<KoShapeReorderCommand::IndexedShape> indexedShapes;
                std::transform(shapes.constBegin(), shapes.constEnd(),
                               std::back_inserter(indexedShapes),
                    [zIndexOffset] (KoShape *shape) {
                        KoShapeReorderCommand::IndexedShape indexedShape(shape);
                        indexedShape.zIndex += zIndexOffset;
                        return indexedShape;
                    });

                indexedShapes = KoShapeReorderCommand::homogenizeZIndexesLazy(indexedShapes);

                KoShapeReorderCommand cmd(indexedShapes);
                cmd.redo();
            }

            KUndo2Command *parentCommand = new KUndo2Command(kundo2_i18n("Paste shapes"));

            const bool forceCreateNewLayer = flags.testFlag(KisPasteActionFactory::ForceNewLayer);
            KoShapeContainer *parentLayer = nullptr;

            if (forceCreateNewLayer) {
                parentLayer = canvas->shapeController()->documentBase()->createParentForShapes(shapes, true, parentCommand);
            }

            canvas->shapeController()->addShapesDirect(shapes, parentLayer, parentCommand);

            QPointF finalShapesOffset;


            if (flags.testFlag(KisPasteActionFactory::PasteAtCursor)) {
                QRectF boundingRect = KoShape::boundingRect(shapes);
                const QPointF cursorPos = canvas->canvasController()->currentCursorPosition();
                finalShapesOffset = cursorPos - boundingRect.center();

            } else if (!forceCreateNewLayer) {
                bool foundOverlapping = false;

                QRectF boundingRect = KoShape::boundingRect(shapes);
                const QPointF offsetStep = 0.1 * QPointF(boundingRect.width(), boundingRect.height());

                QPointF offset;

                Q_FOREACH (KoShape *shape, shapes) {
                    QRectF br1 = shape->boundingRect();

                    bool hasOverlappingShape = false;

                    do {
                        hasOverlappingShape = false;

                        // we cannot use shapesAt() here, because the groups are not
                        // handled in the shape manager's tree
                        QList<KoShape*> conflicts = shapeManager->shapes();

                        Q_FOREACH (KoShape *intersectedShape, conflicts) {
                            if (intersectedShape == shape) continue;

                            QRectF br2 = intersectedShape->boundingRect();

                            const qreal tolerance = 2.0; /* pt */
                            if (KisAlgebra2D::fuzzyCompareRects(br1, br2, tolerance)) {
                                br1.translate(offsetStep.x(), offsetStep.y());
                                offset += offsetStep;

                                hasOverlappingShape = true;
                                foundOverlapping = true;
                                break;
                            }
                        }
                    } while (hasOverlappingShape);

                    if (foundOverlapping) break;
                }

                if (foundOverlapping) {
                    finalShapesOffset = offset;
                }
            }

            const QRectF documentRect = canvas->shapeController()->documentRect();

            if (!forceCreateNewLayer) {
                finalShapesOffset += getFittingOffset(shapes, finalShapesOffset, documentRect, 0.1);
            }

            if (!finalShapesOffset.isNull()) {
                new KoShapeMoveCommand(shapes, finalShapesOffset, parentCommand);
            }

            canvas->addCommand(parentCommand);

            Q_FOREACH (KoShape *shape, shapes) {
                canvas->selectedShapesProxy()->selection()->select(shape);
            }

            result = true;
        }
    }

    return result;
}

}

void KisPasteActionFactory::run(Flags flags, KisViewManager *view)
{
    bool pasteAtCursorPosition = flags.testFlag(PasteAtCursor);

    KisImageSP image = view->image();
    if (!image) return;

    const QPointF docPos = view->canvasBase()->canvasController()->currentCursorPosition();

    // Activate the current tool's paste functionality
    if (view->canvasBase()->toolProxy()->paste()) {
        // XXX: "Add saving of XML data for Paste of shapes"
        return;
    }

    // Paste shapes
    if (tryPasteShapes(flags, view)) {
        return;
    }

    const QRect fittingBounds =
        pasteAtCursorPosition ? QRect() : image->bounds();

    // If no shapes, check for layers
    if (KisClipboard::instance()->hasLayers()) {
        const QPointF offsetTopLeft = [&]() -> QPointF {
            KisPaintDeviceSP clip =
                KisClipboard::instance()->clipFromKritaLayers(
                    image->colorSpace());

            if (!clip) {
                pasteAtCursorPosition = false;
                return {};
            }


            QPointF imagePos;
            if (pasteAtCursorPosition) {
                imagePos =
                    view->canvasBase()->coordinatesConverter()->documentToImage(
                        docPos);

            } else if (!clip->exactBounds().intersects(image->bounds())) {
                 // BUG:459111
                pasteAtCursorPosition = true;
                imagePos = QPointF(image->bounds().center());
            } else {
                return {};
            }
            const QPointF offset =
                (imagePos - QRectF(clip->exactBounds()).center()).toPoint();
            return offset;
        }();


        if (view->selection()) {
            /// TODO: we are relying on a sticky translated string from KisImageLayerAddCommand
            ///       change the string after Krita 5.2.5
            KisProcessingApplicator *ap = beginAction(view, kundo2_i18n("Add Layer"));
            KUndo2Command *deselectCmd = new KisDeselectActiveSelectionCommand(view->selection(), image);
            ap->applyCommand(deselectCmd);
            view->nodeManager()->pasteLayersFromClipboard(pasteAtCursorPosition,
                                                          offsetTopLeft,
                                                          ap);
            endAction(ap, KisOperationConfiguration(id()).toXML());
        } else {
            view->nodeManager()->pasteLayersFromClipboard(pasteAtCursorPosition,
                                                          offsetTopLeft);
        }

        return;
    }

    KisTimeSpan range;
    KisPaintDeviceSP clip = KisClipboard::instance()->clip(fittingBounds, true, -1, &range);

    if (clip) {
        if (pasteAtCursorPosition) {
            const QPointF imagePos = view->canvasBase()->coordinatesConverter()->documentToImage(docPos);

            const QPoint offset =
                (imagePos - QRectF(clip->exactBounds()).center()).toPoint();

            clip->setX(clip->x() + offset.x());
            clip->setY(clip->y() + offset.y());
        }

        KisImportCatcher::adaptClipToImageColorSpace(clip, image);
        bool renamePastedLayers = KisConfig(true).renamePastedLayers();
        QString pastedLayerName = renamePastedLayers ? image->nextLayerName() + " " + i18n("(pasted)") :
                                                       image->nextLayerName();
        KisPaintLayerSP newLayer = new KisPaintLayer(image.data(),
                                                     pastedLayerName,
                                                     OPACITY_OPAQUE_U8);
        KisNodeSP aboveNode = view->activeLayer();
        KisNodeSP parentNode = aboveNode ? aboveNode->parent() : image->root();

        if (range.isValid()) {
            newLayer->enableAnimation();
            KisKeyframeChannel *channel = newLayer->getKeyframeChannel(KisKeyframeChannel::Raster.id(), true);
            KisRasterKeyframeChannel *rasterChannel = dynamic_cast<KisRasterKeyframeChannel*>(channel);
            KIS_SAFE_ASSERT_RECOVER_RETURN(rasterChannel);
            rasterChannel->importFrame(range.start(), clip, nullptr);

            if (!range.isInfinite()) {
                rasterChannel->addKeyframe(range.end() + 1, 0);
            }
        } else {
            const QRect rc = clip->extent();
            KisPainter::copyAreaOptimized(rc.topLeft(), clip, newLayer->paintDevice(), rc);
        }

        KUndo2Command *cmd = new KisImageLayerAddCommand(image, newLayer, parentNode, aboveNode);
        KisProcessingApplicator *ap = beginAction(view, cmd->text());
        ap->applyCommand(cmd, KisStrokeJobData::SEQUENTIAL, KisStrokeJobData::NORMAL);
        
        if (view->selection()) {
            KUndo2Command *deselectCmd = new KisDeselectActiveSelectionCommand(view->selection(), image);
            ap->applyCommand(deselectCmd, KisStrokeJobData::SEQUENTIAL, KisStrokeJobData::NORMAL);
        }

        if (KisConfig(true).activateTransformToolAfterPaste()) {
            KUndo2Command *transformToolCmd = new KisTransformToolActivationCommand(view);
            ap->applyCommand(transformToolCmd, KisStrokeJobData::BARRIER, KisStrokeJobData::NORMAL);
        }
        
        endAction(ap, KisOperationConfiguration(id()).toXML());
    }
}

void KisPasteIntoActionFactory::run(KisViewManager *viewManager)
{
    if (!viewManager->activeDevice()) return;

    KisImageSP image = viewManager->image();
    if (!image) return;

    QRect imageBounds = image->bounds();

    KisPaintDeviceSP clipdev = KisClipboard::instance()->clipFromKritaLayers(image->colorSpace());
    KisPaintDeviceSP clip = clipdev ? new KisPaintDevice(*clipdev) : nullptr;

    if (clip)
    {
        QRect clipBounds = clip->exactBounds();

        if (!clipBounds.intersects(imageBounds))
        {
            QPoint diff = imageBounds.center() - clipBounds.center();
            clip->setX(diff.x());
            clip->setY(diff.y());
        }
    }
    else
    {
        clip = KisClipboard::instance()->clip(imageBounds, true, -1, nullptr);
    }

    if (!clip) return;

    KisImportCatcher::adaptClipToImageColorSpace(clip, image);

    if (viewManager->selection()) {
        KUndo2Command *deselectCmd = new KisDeselectActiveSelectionCommand(viewManager->selection(), image);
        KisProcessingApplicator::runSingleCommandStroke(viewManager->image(), deselectCmd);
    }

    KisTool* tool = dynamic_cast<KisTool*>(KoToolManager::instance()->toolById(viewManager->canvasBase(), "KisToolTransform"));
    KIS_ASSERT(tool);
    tool->newActivationWithExternalSource(clip);
}

void KisPasteNewActionFactory::run(KisViewManager *viewManager)
{
    Q_UNUSED(viewManager);

    KisPaintDeviceSP clip = KisClipboard::instance()->clip(QRect(), true);
    if (!clip) return;

    QRect rect = clip->exactBounds();
    if (rect.isEmpty()) return;

    KisDocument *doc = KisPart::instance()->createDocument();
    doc->documentInfo()->setAboutInfo("title", i18n("Untitled"));
    KisImageSP image = new KisImage(doc->createUndoStore(),
                                    rect.width(),
                                    rect.height(),
                                    clip->colorSpace(),
                                    i18n("Pasted"));
    bool renamePastedLayers = KisConfig(true).renamePastedLayers();
    QString pastedLayerName = renamePastedLayers ? image->nextLayerName() + " " + i18n("(pasted)") :
                                                   image->nextLayerName();
    KisPaintLayerSP layer =
            new KisPaintLayer(image.data(), pastedLayerName,
                              OPACITY_OPAQUE_U8, clip->colorSpace());

    KisPainter::copyAreaOptimized(QPoint(), clip, layer->paintDevice(), rect);

    image->addNode(layer.data(), image->rootLayer());
    doc->setCurrentImage(image);
    KisPart::instance()->addDocument(doc);

    KisMainWindow *win = viewManager->mainWindow();
    win->addViewAndNotifyLoadingCompleted(doc);
}

void KisPasteReferenceActionFactory::run(KisViewManager *viewManager)
{
    KisCanvas2 *canvasBase = viewManager->canvasBase();
    if (!canvasBase) return;

    KisReferenceImage* reference = KisReferenceImage::fromClipboard(*canvasBase->coordinatesConverter());
    if (!reference) return;

    KisDocument *doc = viewManager->document();
    canvasBase->addCommand(KisReferenceImagesLayer::addReferenceImages(doc, {reference}));

    KoToolManager::instance()->switchToolRequested("ToolReferenceImages");
}

void KisPasteShapeStyleActionFactory::run(KisViewManager *view)
{
    KoSvgPaste paste;

    KisCanvas2 *canvas = view->canvasBase();

    KoShapeManager *shapeManager = canvas->shapeManager();
    QList<KoShape*> selectedShapes = shapeManager->selection()->selectedEditableShapes();

    if (selectedShapes.isEmpty()) return;

    if (paste.hasShapes()) {
        KoCanvasBase *canvas = view->canvasBase();

        QSizeF fragmentSize;
        QList<KoShape*> shapes =
            paste.fetchShapes(canvas->shapeController()->documentRectInPixels(),
                              canvas->shapeController()->pixelsPerInch(), &fragmentSize);

        if (!shapes.isEmpty()) {
            KoShape *referenceShape = shapes.first();


            KUndo2Command *parentCommand = new KUndo2Command(kundo2_i18n("Paste Style"));

            new KoShapeBackgroundCommand(selectedShapes, referenceShape->background(), parentCommand);
            new KoShapeStrokeCommand(selectedShapes, referenceShape->stroke(), parentCommand);


            canvas->addCommand(parentCommand);
        }

        qDeleteAll(shapes);
    }
}
