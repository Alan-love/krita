/*
 *  SPDX-FileCopyrightText: 2011 Dmitry Kazakov <dimula73@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "kis_painting_information_builder.h"

#include <KoPointerEvent.h>

#include "kis_config.h"
#include "kis_config_notifier.h"

#include "kis_cubic_curve.h"
#include "kis_speed_smoother.h"

#include <KoCanvasResourceProvider.h>
#include "kis_canvas_resource_provider.h"


/***********************************************************************/
/*           KisPaintingInformationBuilder                             */
/***********************************************************************/


const int KisPaintingInformationBuilder::LEVEL_OF_PRESSURE_RESOLUTION = 1024;


KisPaintingInformationBuilder::KisPaintingInformationBuilder()
    : m_speedSmoother(new KisSpeedSmoother()),
      m_pressureDisabled(false)
{
    connect(KisConfigNotifier::instance(), SIGNAL(configChanged()),
            SLOT(updateSettings()));

    updateSettings();
}

KisPaintingInformationBuilder::~KisPaintingInformationBuilder()
{

}

void KisPaintingInformationBuilder::updateSettings()
{
    KisConfig cfg(true);
    const KisCubicCurve curve(cfg.pressureTabletCurve());
    m_pressureSamples = curve.floatTransfer(LEVEL_OF_PRESSURE_RESOLUTION + 1);
    m_maxAllowedSpeedValue = cfg.readEntry("maxAllowedSpeedValue", 30);

    // The setting is stored in [-180, 180] degrees range to match the UI.
    // We now need to convert it to [0, 360) range.
    m_tiltDirectionOffset = cfg.readEntry("tiltDirectionOffset", 0);
    if (m_tiltDirectionOffset < 0.0) {
        m_tiltDirectionOffset += 360.0;
    }

    m_speedSmoother->updateSettings();
}

KisPaintInformation KisPaintingInformationBuilder::startStroke(KoPointerEvent *event,
                                                               int timeElapsed,
                                                               const KoCanvasResourceProvider *manager)
{
    if (manager) {
        m_pressureDisabled = manager->resource(KoCanvasResource::DisablePressure).toBool();
    }

    m_startPoint = event->point;
    return createPaintingInformation(event, timeElapsed);

}

KisPaintInformation KisPaintingInformationBuilder::continueStroke(KoPointerEvent *event,
                                                                  int timeElapsed)
{
    return createPaintingInformation(event, timeElapsed);
}

QPointF KisPaintingInformationBuilder::adjustDocumentPoint(const QPointF &point, const QPointF &/*startPoint*/)
{
    return point;
}

QPointF KisPaintingInformationBuilder::documentToImage(const QPointF &point)
{
    return point;
}

QPointF KisPaintingInformationBuilder::imageToDocument(const QPointF &point)
{
    return point;
}

QPointF KisPaintingInformationBuilder::imageToView(const QPointF &point)
{
    return point;
}

qreal KisPaintingInformationBuilder::calculatePerspective(const QPointF &documentPoint)
{
    Q_UNUSED(documentPoint);
    return 1.0;
}

qreal KisPaintingInformationBuilder::canvasRotation() const
{
    return 0;
}

bool KisPaintingInformationBuilder::canvasMirroredX() const
{
    return false;
}

bool KisPaintingInformationBuilder::canvasMirroredY() const
{
    return false;
}

KisPaintInformation KisPaintingInformationBuilder::createPaintingInformation(KoPointerEvent *event,
                                                                             int timeElapsed)
{

    QPointF adjusted = adjustDocumentPoint(event->point, m_startPoint);
    QPointF imagePoint = documentToImage(adjusted);
    qreal perspective = calculatePerspective(adjusted);
    const qreal speed = m_speedSmoother->getNextSpeed(imageToView(imagePoint), event->time());

    KisPaintInformation pi(imagePoint,
                           !m_pressureDisabled ? 1.0 : pressureToCurve(event->pressure()),
                           event->xTilt(), event->yTilt(),
                           event->rotation(),
                           event->tangentialPressure(),
                           perspective,
                           timeElapsed,
                           qMin(1.0, speed / qreal(m_maxAllowedSpeedValue)));

    pi.setCanvasRotation(canvasRotation());
    pi.setCanvasMirroredH(canvasMirroredX());
    pi.setCanvasMirroredV(canvasMirroredY());
    pi.setTiltDirectionOffset(m_tiltDirectionOffset);

    return pi;
}

KisPaintInformation KisPaintingInformationBuilder::hover(const QPointF &imagePoint,
                                                         const KoPointerEvent *event,
                                                         bool isStrokeStarted)
{
    const qreal perspective = calculatePerspective(imageToDocument(imagePoint));

    const qreal speed = !isStrokeStarted && event ?
        m_speedSmoother->getNextSpeed(imageToView(imagePoint), event->time()) :
        m_speedSmoother->lastSpeed();

    if (event) {
        return KisPaintInformation::createHoveringModeInfo(imagePoint,
                                                           PRESSURE_DEFAULT,
                                                           event->xTilt(), event->yTilt(),
                                                           event->rotation(),
                                                           event->tangentialPressure(),
                                                           perspective,
                                                           qMin(1.0, speed / qreal(m_maxAllowedSpeedValue)),
                                                           canvasRotation(),
                                                           canvasMirroredX(),
                                                           canvasMirroredY(),
                                                           m_tiltDirectionOffset);
    } else {
        KisPaintInformation pi = KisPaintInformation::createHoveringModeInfo(imagePoint);
        pi.setCanvasRotation(canvasRotation());
        pi.setCanvasMirroredH(canvasMirroredX());
        pi.setCanvasMirroredV(canvasMirroredY());
        pi.setTiltDirectionOffset(m_tiltDirectionOffset);
        return pi;
    }
}

qreal KisPaintingInformationBuilder::pressureToCurve(qreal pressure)
{
    return KisCubicCurve::interpolateLinear(pressure, m_pressureSamples);
}

void KisPaintingInformationBuilder::reset()
{
    m_speedSmoother->clear();
}

/***********************************************************************/
/*           KisConverterPaintingInformationBuilder                        */
/***********************************************************************/

#include "kis_coordinates_converter.h"

KisConverterPaintingInformationBuilder::KisConverterPaintingInformationBuilder(const KisCoordinatesConverter *converter)
    : m_converter(converter)
{
}

QPointF KisConverterPaintingInformationBuilder::documentToImage(const QPointF &point)
{
    return m_converter->documentToImage(point);
}

QPointF KisConverterPaintingInformationBuilder::imageToDocument(const QPointF &point)
{
    return m_converter->imageToDocument(point);
}

QPointF KisConverterPaintingInformationBuilder::imageToView(const QPointF &point)
{
    return m_converter->imageToWidget(point);
}

qreal KisConverterPaintingInformationBuilder::canvasRotation() const
{
    return m_converter->rotationAngle();
}

bool KisConverterPaintingInformationBuilder::canvasMirroredX() const
{
    return m_converter->xAxisMirrored();
}

bool KisConverterPaintingInformationBuilder::canvasMirroredY() const
{
    return m_converter->yAxisMirrored();
}

/***********************************************************************/
/*           KisToolFreehandPaintingInformationBuilder                 */
/***********************************************************************/

#include "kis_tool_freehand.h"
#include "kis_canvas2.h"

KisToolFreehandPaintingInformationBuilder::KisToolFreehandPaintingInformationBuilder(KisToolFreehand *tool)
    : m_tool(tool)
{
}

QPointF KisToolFreehandPaintingInformationBuilder::documentToImage(const QPointF &point)
{
    return m_tool->convertToPixelCoord(point);
}

QPointF KisToolFreehandPaintingInformationBuilder::imageToDocument(const QPointF &point)
{
    KisCanvas2 *canvas = dynamic_cast<KisCanvas2*>(m_tool->canvas());
    KIS_ASSERT_RECOVER_RETURN_VALUE(canvas, point);
    return canvas->coordinatesConverter()->imageToDocument(point);
}

QPointF KisToolFreehandPaintingInformationBuilder::imageToView(const QPointF &point)
{
    return m_tool->pixelToView(point);
}

QPointF KisToolFreehandPaintingInformationBuilder::adjustDocumentPoint(const QPointF &point, const QPointF &startPoint)
{
    return m_tool->adjustPosition(point, startPoint);
}

qreal KisToolFreehandPaintingInformationBuilder::calculatePerspective(const QPointF &documentPoint)
{
    return m_tool->calculatePerspective(documentPoint);
}

qreal KisToolFreehandPaintingInformationBuilder::canvasRotation() const
{
    KisCanvas2 *canvas = dynamic_cast<KisCanvas2*>(m_tool->canvas());
    KIS_ASSERT_RECOVER_RETURN_VALUE(canvas, 0.0);
    return canvas->coordinatesConverter()->rotationAngle();
}

bool KisToolFreehandPaintingInformationBuilder::canvasMirroredX() const
{
    KisCanvas2 *canvas = dynamic_cast<KisCanvas2*>(m_tool->canvas());
    KIS_ASSERT_RECOVER_RETURN_VALUE(canvas, false);
    return canvas->coordinatesConverter()->xAxisMirrored();
}

bool KisToolFreehandPaintingInformationBuilder::canvasMirroredY() const
{
    KisCanvas2 *canvas = dynamic_cast<KisCanvas2*>(m_tool->canvas());
    KIS_ASSERT_RECOVER_RETURN_VALUE(canvas, false);
    return canvas->coordinatesConverter()->yAxisMirrored();
}
