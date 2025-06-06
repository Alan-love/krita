/* This file is part of the Krita project
 *
 * SPDX-FileCopyrightText: 2014 Boudewijn Rempt <boud@valdyas.org>
 *
 * SPDX-License-Identifier: LGPL-2.0-or-later
 */

#ifndef KISDOCUMENT_H
#define KISDOCUMENT_H

#include <QDateTime>
#include <QList>
#include <QFileInfo>

#include <klocalizedstring.h>

#include <kundo2stack.h>
#include <KoColorSet.h>

#include <kis_image.h>
#include <KisImportExportFilter.h>
#include <kis_properties_configuration.h>
#include <kis_types.h>
#include <kis_painting_assistant.h>
#include <KisReferenceImage.h>
#include <kis_debug.h>
#include <KisImportExportUtils.h>
#include <kis_config.h>
#include "StoryboardItem.h"

#include "kritaui_export.h"

#include <memory>

class QString;

class KUndo2Command;
class KoUnit;

class KoColor;
class KoColorSpace;
class KoShapeControllerBase;
class KoShapeLayer;
class KoStore;
class KoDocumentInfo;
class KoDocumentInfoDlg;
class KisImportExportManager;
class KisUndoStore;
class KisPart;
class KisGridConfig;
class KisGuidesConfig;
class KisMirrorAxisConfig;
class QDomDocument;
class KisReferenceImagesLayer;

#define KIS_MIME_TYPE "application/x-krita"

/**
 *  KisDocument contains the image and keeps track of loading,
 *  modification, undo stack and saving.
 */
class KRITAUI_EXPORT KisDocument : public QObject
{
    Q_OBJECT
protected:

    explicit KisDocument(bool addStorage = true);

    /**
     * @brief KisDocument makes a deep copy of the document \p rhs.
     *        The caller *must* ensure that the image is properly
     *        locked and is in consistent state before asking for
     *        cloning.
     * @param rhs the source document to copy from
     */
    explicit KisDocument(const KisDocument &rhs, bool addStorage);

public:
    enum OpenFlag {
        None = 0,
        DontAddToRecent = 0x1,
        RecoveryFile = 0x2
    };
    Q_DECLARE_FLAGS(OpenFlags, OpenFlag)

    /**
     *  Destructor.
     *
     * The destructor does not delete any attached KisView objects and it does not
     * delete the attached widget as returned by widget().
     */
    ~KisDocument();

    /**
     * Temporary storage for the resources that are embedded into other
     * resources used by the document. E.g. patterns embedded into layer
     * styles.
     */
    QString embeddedResourcesStorageId() const;

    /**
     * Persistent storage for the resources that are linked but the resources
     * embedded in the document. These resources are not embedded into their own
     * container resource, so they should be stored by the document
     *
     * All these resources are saved into the document itself and loaded
     * alongside the document.
     */
    QString linkedResourcesStorageId() const;

    /**
     * @brief creates a clone of the document and returns it. Please make sure that you
     * hold all the necessary locks on the image before asking for a clone!
     */
    KisDocument *clone(bool addStorage = false);

    /**
     * @brief openPath Open a Path
     * @param path Path to file
     * @param flags Control specific behavior
     * @return success status
     */
    bool openPath(const QString &path, OpenFlags flags = None);

    /**
     * Opens the document given by @p path, without storing the Path
     * in the KisDocument.
     * Call this instead of openPath() to implement KisMainWindow's
     * File --> Import feature.
     *
     * @note This will call openPath(). To differentiate this from an ordinary
     *       Open operation (in any reimplementation of openPath() or openFile())
     *       call isImporting().
     */
    bool importDocument(const QString &path);

    /**
     * Saves the document as @p path without changing the state of the
     * KisDocument (Path, modified flag etc.). Call this instead of
     * saveAs() to implement KisMainWindow's File --> Export feature.
     * Make sure to provide two separate bool parameters otherwise it will mix them
     */
    bool exportDocument(const QString &path, const QByteArray &mimeType,bool isAdvancedExporting = false, bool showWarnings = false, KisPropertiesConfigurationSP exportConfiguration = 0);
    /**
     * Exports he document is a synchronous way. The caller must ensure that the
     * image is not accessed by any other actors, because the exporting happens
     * without holding the image lock.
     */
    bool exportDocumentSync(const QString &path, const QByteArray &mimeType, KisPropertiesConfigurationSP exportConfiguration = 0);

private:
    bool exportDocumentImpl(const KritaUtils::ExportFileJob &job, KisPropertiesConfigurationSP exportConfiguration, bool isAdvancedExporting= false);

public:
    /**
     * @brief Sets whether the document can be edited or is read only.
     *
     * This recursively applied to all child documents and
     * KisView::updateReadWrite is called for every attached
     * view.
     */
    void setReadWrite(bool readwrite = true);

    /**
     * To be preferred when a document exists. It is fast when calling
     * it multiple times since it caches the result that readNativeFormatMimeType()
     * delivers.
     * This comes from the X-KDE-NativeMimeType key in the .desktop file.
     */
    static QByteArray nativeFormatMimeType() { return KIS_MIME_TYPE; }

    /// Checks whether a given mimeType can be handled natively.
    bool isNativeFormat(const QByteArray& mimeType) const;

    /// Returns a list of the mimeTypes considered "native", i.e. which can
    /// be saved by KisDocument without a filter, in *addition* to the main one
    static QStringList extraNativeMimeTypes() { return QStringList() << KIS_MIME_TYPE; }

    /**
     * Returns the actual mimeType of the document
     */
    QByteArray mimeType() const;

    /**
     * @brief Sets the mime type for the document.
     *
     * When choosing "save as" this is also the mime type
     * selected by default.
     */
    void setMimeType(const QByteArray & mimeType);

    /**
     * @return true if file operations should inhibit the option dialog
     */
    bool fileBatchMode() const;

    /**
     * @param batchMode if true, do not show the option dialog for file operations.
     */
    void setFileBatchMode(const bool batchMode);

    /**
     * Sets the error message to be shown to the user (use i18n()!)
     * when loading or saving fails.
     * If you asked the user about something and they chose "Cancel",
     */
    void setErrorMessage(const QString& errMsg);

    /**
     * Return the last error message. Usually KisDocument takes care of
     * showing it; this method is mostly provided for non-interactive use.
     */
    QString errorMessage() const;

    /**
     * Sets the warning message to be shown to the user (use i18n()!)
     * when loading or saving fails.
     */
    void setWarningMessage(const QString& warningMsg);

    /**
     * Return the last warning message set by loading or saving. Warnings
     * mean that the document could not be completely loaded, but the errors
     * were not absolutely fatal.
     */
    QString warningMessage() const;

    /**
     * @brief Generates a preview picture of the document
     * @note The preview is used in the File Dialog and also to create the Thumbnail
     */
    QPixmap generatePreview(const QSize& size);

    /**
     *  @brief Sets the document to empty.
     *
     *  Used after loading a template
     *  (which is not empty, but not the user's input).
     *
     *  @see isEmpty()
     */
    void setEmpty(bool empty = true);

    /**
     *  Return a correctly created QDomDocument for this KisDocument,
     *  including processing instruction, complete DOCTYPE tag (with systemId and publicId), and root element.
     *  @param tagName the name of the tag for the root element
     *  @param version the DTD version (usually the application's version).
     */
    QDomDocument createDomDocument(const QString& tagName, const QString& version) const;

    /**
     *  Return a correctly created QDomDocument for an old (1.3-style) Krita document,
     *  including processing instruction, complete DOCTYPE tag (with systemId and publicId), and root element.
     *  This static method can be used e.g. by filters.
     *  @param appName the app's instance name, e.g. words, kspread, kpresenter etc.
     *  @param tagName the name of the tag for the root element, e.g. DOC for words/kpresenter.
     *  @param version the DTD version (usually the application's version).
     */
    static QDomDocument createDomDocument(const QString& appName, const QString& tagName, const QString& version);

   /**
     *  Loads a document in the native format from a given Path.
     *  Reimplement if your native format isn't XML.
     *
     *  @param file the file to load - usually KReadOnlyPart::m_file or the result of a filter
     */
    bool loadNativeFormat(const QString & file);

    /**
     * allow to activate or deactivate autosave on document, independently of auto save delay
     *
     * the value is independant of auto save delay
     */
    void setAutoSaveActive(bool autoSaveIsActive);

    /**
     * indicate if autosave is active or inactive
     *
     * the value is independant of auto save delay
     */
    bool isAutoSaveActive();

    /**
     * Set standard autosave interval that is set by a config file
     */
    void setNormalAutoSaveInterval();

    /**
     * Set emergency interval that autosave uses when the image is busy,
     * by default it is 10 sec
     */
    void setEmergencyAutoSaveInterval();

    /**
     * Disable autosave
     */
    void setInfiniteAutoSaveInterval();

    /**
     * @return the information concerning this document.
     * @see KoDocumentInfo
     */
    KoDocumentInfo *documentInfo() const;

    /**
     * Performs a cleanup of unneeded backup files
     */
    void removeAutoSaveFiles(const QString &autosaveBaseName, bool wasRecovered);

    /**
     * Returns true if this document or any of its internal child documents are modified.
     */
    bool isModified() const;

    /**
     * @return caption of the document
     *
     * Caption is of the form "[title] - [path]",
     * built out of the document info (title) and pretty-printed
     * document Path.
     * If the title is not present, only the Path it returned.
     */
    QString caption() const;

    /**
     * Sets the document Path to empty Path
     * After using loadNativeFormat on a template, one wants
     * to set the path to QString()
     */
    void resetPath();

    /**
     * @internal (public for KisMainWindow)
     */
    void setMimeTypeAfterLoading(const QString& mimeType);

    /**
     * Returns the unit used to display all measures/distances.
     */
    KoUnit unit() const;

    /**
     * Sets the unit used to display all measures/distances.
     */
    void setUnit(const KoUnit &unit);

    KisGridConfig gridConfig() const;
    void setGridConfig(const KisGridConfig &config);

    /// returns the guides data for this document.
    const KisGuidesConfig& guidesConfig() const;
    void setGuidesConfig(const KisGuidesConfig &data);

    /**
     * @brief linkedDocumentResources List returns all the resources
     * linked to the document, such as palettes
     *
     * In some cases (e.g. when the document is temporary), the
     * underlying document storage will not be registered in the
     * resource system, so we cannot get fully initialized resources
     * from it (resourceId(), active(), md5() and storageLocation()
     * fields will be uninitialized). Therefore we just return
     * KoEmbeddedResource which is suitable for saving this data into
     * hard drive.
     *
     * The returned KoResourceLoadResult object can either be in
     * EmbeddedResource or FailedLink state. The former means the
     * resource has been prepared for embedding, the latter means
     * there was some issue with serializing the resource. In the
     * latter case the called should check result.signature() to
     * find out which resource has failed.
     *
     * NOTE: the returned result can **NOT** have ExistingResource state!
     */
    QList<KoResourceLoadResult> linkedDocumentResources();

    /**
     * @brief setPaletteList replaces the palettes in the document's local resource storage with the list
     * of palettes passed to this function. It will then Q_EMIT sigPaletteListChanged with both the old and
     * the new list, if emitSignal is true.
     */
    void setPaletteList(const QList<KoColorSetSP> &paletteList, bool emitSignal = false);

    /**
     * @brief returns the list of pointers to storyboard Items for the document
     */
    StoryboardItemList getStoryboardItemList();

    /**
     * @brief sets the storyboardItemList in the document, emits empty signal if emitSignal is true.
     */
    void setStoryboardItemList(const StoryboardItemList &storyboardItemList, bool emitSignal = false);

    /**
     * @brief returns the list of comments for the storyboard docker in the document
     */
    QVector<StoryboardComment> getStoryboardCommentsList();

    /**
     * @brief sets the  list of comments for the storyboard docker in the document, emits empty signal if emitSignal is true.
     */
    void setStoryboardCommentList(const QVector<StoryboardComment> &storyboardCommentList, bool emitSignal = false);

    QVector<QFileInfo> getAudioTracks() const;
    void setAudioTracks(QVector<QFileInfo> f);

    void setAudioVolume(qreal level);
    qreal getAudioLevel();

    const KisMirrorAxisConfig& mirrorAxisConfig() const;
    void setMirrorAxisConfig(const KisMirrorAxisConfig& config);

    void clearUndoHistory();

    /**
     *  Sets the modified flag on the document. This means that it has
     *  to be saved or not before deleting it.
     */
    void setModified(bool _mod);

    void setRecovered(bool value);
    bool isRecovered() const;

    void updateEditingTime(bool forceStoreElapsed);

    /**
     * Returns the global undo stack
     */
    KUndo2Stack *undoStack();


    /**
     * @brief importExportManager gives access to the internal import/export manager
     * @return the document's import/export manager
     */
    KisImportExportManager *importExportManager() const;

    /**
     * @brief serializeToNativeByteArray daves the document into a .kra file written
     * to a memory-based byte-array
     * @return a byte array containing the .kra file
     */
    QByteArray serializeToNativeByteArray();


    /**
     * @brief isInSaving shown if the document has any (background) saving process or not
     * @return true if there is some saving in action
     */
    bool isInSaving() const;

Q_SIGNALS:

    /**
     * This signal is emitted when the unit is changed by setUnit().
     * It is common to connect views to it, in order to change the displayed units
     * (e.g. in the rulers)
     */
    void unitChanged(const KoUnit &unit);

    /**
     * Emitted e.g. at the beginning of a save operation
     * This is emitted by KisDocument and used by KisView to display a statusbar message
     */
    void statusBarMessage(const QString& text, int timeout = 0);

    /**
     * Emitted e.g. at the end of a save operation
     * This is emitted by KisDocument and used by KisView to clear the statusbar message
     */
    void clearStatusBarMessage();

    /**
    * Emitted when the document is modified
    */
    void modified(bool);

    void sigReadWriteChanged(bool value);
    void sigRecoveredChanged(bool value);
    void sigPathChanged(const QString &path);

    void sigLoadingFinished();

    void sigSavingFinished(const QString &filePath);

    void sigGuidesConfigChanged(const KisGuidesConfig &config);

    void sigBackgroundSavingFinished(KisImportExportErrorCode status, const QString &errorMessage, const QString &warningMessage);

    void sigCompleteBackgroundSaving(const KritaUtils::ExportFileJob &job, KisImportExportErrorCode status, const QString &errorMessage, const QString &warningMessage);

    void sigReferenceImagesChanged();

    void sigMirrorAxisConfigChanged();

    void sigGridConfigChanged(const KisGridConfig &config);

    void sigReferenceImagesLayerChanged(KisSharedPtr<KisReferenceImagesLayer> layer);

    /**
     * Emitted when the palette list has changed.
     * The pointers in oldPaletteList are to be deleted by the resource server.
     **/
    void sigPaletteListChanged(const QList<KoColorSetSP> &oldPaletteList, const QList<KoColorSetSP> &newPaletteList);

    void sigAssistantsChanged();

    void sigStoryboardItemListChanged();

    void sigStoryboardCommentListChanged();

    void sigAudioTracksChanged();

    void sigAudioLevelChanged(qreal level);

private Q_SLOTS:
    void finishExportInBackground();
    void slotChildCompletedSavingInBackground(KisImportExportErrorCode status, const QString &errorMessage, const QString &warningMessage);
    void slotCompleteAutoSaving(const KritaUtils::ExportFileJob &job, KisImportExportErrorCode status, const QString &errorMessage, const QString &warningMessage);

    void slotCompleteSavingDocument(const KritaUtils::ExportFileJob &job, KisImportExportErrorCode status, const QString &errorMessage, const QString &warningMessage);

    void slotInitiateAsyncAutosaving(KisDocument *clonedDocument);
    void slotDocumentCloningCancelled();

    void slotPerformIdleRoutines();

private:

    friend class KisPart;
    friend class SafeSavingLocker;

    KritaUtils::BackgroudSavingStartResult initiateSavingInBackground(const QString actionName,
                                    const QObject *receiverObject, const char *receiverMethod,
                                    const KritaUtils::ExportFileJob &job,
                                    KisPropertiesConfigurationSP exportConfiguration,
                                    std::unique_ptr<KisDocument> &&optionalClonedDocument, bool isAdvancedExporting = false);

    KritaUtils::BackgroudSavingStartResult initiateSavingInBackground(const QString actionName,
                                    const QObject *receiverObject, const char *receiverMethod,
                                    const KritaUtils::ExportFileJob &job,
                                    KisPropertiesConfigurationSP exportConfiguration, bool isAdvancedExporting =false );

    bool startExportInBackground(const QString &actionName, const QString &location,
                                 const QString &realLocation,
                                 const QByteArray &mimeType,
                                 bool showWarnings,
                                 KisPropertiesConfigurationSP exportConfiguration, bool isAdvancedExporting= false);

    /**
     * Activate/deactivate/configure the autosave feature.
     * @param delay in seconds, 0 to disable
     */
    void setAutoSaveDelay(int delay);

    /**
     * Generate a name for the document.
     */
    QString newObjectName();

    QString generateAutoSaveFileName(const QString & path) const;

    /**
     *  Loads a document
     *
     *  Applies a filter if necessary, and calls loadNativeFormat in any case
     *  You should not have to reimplement, except for very special cases.
     *
     * NOTE: this method also creates a new KisView instance!
     *
     * This method is called from the KReadOnlyPart::openPath method.
     */
    bool openFile();

public:

    bool isAutosaving() const;

public:

    QString localFilePath() const;
    void setLocalFilePath( const QString &localFilePath );

    KoDocumentInfoDlg* createDocumentInfoDialog(QWidget *parent, KoDocumentInfo *docInfo) const;

    bool isReadWrite() const;

    QString path() const;
    void setPath(const QString &path);

    bool closePath(bool promptToSave = true);

    bool saveAs(const QString &path, const QByteArray &mimeType, bool showWarnings, KisPropertiesConfigurationSP exportConfiguration = 0);

    /**
     * Create a new image that has this document as a parent and
     * replace the current image with this image.
     */
    bool newImage(const QString& name, qint32 width, qint32 height, const KoColorSpace * cs, const KoColor &bgColor, KisConfig::BackgroundStyle bgStyle,
                  int numberOfLayers, const QString &imageDescription, const double imageResolution);

    bool isSaving() const;
    void waitForSavingToComplete();


    KisImageWSP image() const;

    /**
     * @brief savingImage provides a detached, shallow copy of the original image that must be used when saving.
     * Any strokes in progress will not be applied to this image, so the result might be missing some data. On
     * the other hand, it won't block.
     *
     * @return a shallow copy of the original image, or 0 is saving is not in progress
     */
    KisImageSP savingImage() const;

    /**
     * Set the current image to the specified image and turn undo on.
     */
    void setCurrentImage(KisImageSP image, bool forceInitialUpdate = true, KisNodeSP preActivatedNode = nullptr);

    /**
     * Set the image of the document preliminary, before the document
     * has completed loading. Some of the document items (shapes) may want
     * to access image properties (bounds and resolution), so we should provide
     * it to them even before the entire image is loaded.
     *
     * Right now, the only use by KoShapeRegistry::createShapeFromOdf(), remove
     * after it is deprecated.
     */
    void hackPreliminarySetImage(KisImageSP image);

    KisUndoStore* createUndoStore();

    /**
     * The shape controller matches internal krita image layers with
     * the flake shape hierarchy.
     */
    KoShapeControllerBase *shapeController() const;

    KoShapeLayer *shapeForNode(KisNodeSP layer) const;

    /**
     * Set the list of nodes that was marked as currently active. Used *only*
     * for saving loading. Never use it for tools or processing.
     */
    void setPreActivatedNode(KisNodeSP activatedNode);

    /**
     * @return the node that was set as active during loading. Used *only*
     * for saving loading. Never use it for tools or processing.
     */
    KisNodeSP preActivatedNode() const;

    /// @return the list of assistants associated with this document
    QList<KisPaintingAssistantSP> assistants() const;

    /// @replace the current list of assistants with @param value
    void setAssistants(const QList<KisPaintingAssistantSP> &value);


    void setAssistantsGlobalColor(QColor color);
    QColor assistantsGlobalColor();

    // Color history if per document (configuration dependent)
    void setColorHistory(const QList<KoColor> &colors);
    QList<KoColor> colorHistory();

    /**
     * Get existing reference images layer or null if none exists.
     */
    KisSharedPtr<KisReferenceImagesLayer> referenceImagesLayer() const;

    void setReferenceImagesLayer(KisSharedPtr<KisReferenceImagesLayer> layer, bool updateImage);

    bool save(bool showWarnings, KisPropertiesConfigurationSP exportConfiguration);

    /**
     * Return the bounding box of the image and associated elements (e.g. reference images)
     */
    QRectF documentBounds() const;

    /**
     * @brief Start saving when android activity is pushed to the background
     */
    void autoSaveOnPause();

Q_SIGNALS:

    void completed();
    void canceled(const QString &);

private Q_SLOTS:

    void setImageModified();
    void setImageModifiedWithoutUndo();

    void slotAutoSave();

    void slotUndoStackCleanChanged(bool value);

    void slotConfigChanged();

    void slotImageRootChanged();


public:
    /**
     * @brief try to clone the image. This method handles all the locking for you. If locking
     *        has failed, no cloning happens
     * @return cloned document on success, null otherwise
     */
    KisDocument *lockAndCloneForSaving();

    KisDocument *lockAndCreateSnapshot();

    void copyFromDocument(const KisDocument &rhs);

private:

    enum CopyPolicy {
        CONSTRUCT = 0, ///< we are copy-constructing a new KisDocument
        REPLACE ///< we are replacing the current KisDocument with another
    };

    void copyFromDocumentImpl(const KisDocument &rhs, CopyPolicy policy);

    QString exportErrorToUserMessage(KisImportExportErrorCode status, const QString &errorMessage);

    QString prettyPath() const;

    bool openPathInternal(const QString &path);

    void slotAutoSaveImpl(std::unique_ptr<KisDocument> &&optionalClonedDocument);

    /// Checks whether we are saving a resource we've been editing, and if so,
    /// uses the resource server to update the resource.
    /// @return true if this was a resource, false if the document needs to be saved
    bool resourceSavingFilter(const QString &path, const QByteArray &mimeType, KisPropertiesConfigurationSP exportConfiguration);

    class Private;
    Private *const d;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(KisDocument::OpenFlags)
Q_DECLARE_METATYPE(KisDocument*)

#endif
