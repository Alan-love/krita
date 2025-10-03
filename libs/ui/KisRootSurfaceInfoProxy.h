/*
 *  SPDX-FileCopyrightText: 2025 Dmitry Kazakov <dimula73@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */
#ifndef KISROOTSURFACEINFOPROXY_H
#define KISROOTSURFACEINFOPROXY_H

#include <kritaui_export.h>

#include <QObject>
#include <QPointer>

class KisSRGBSurfaceColorSpaceManager;
class QWidget;
class QWindow;
class QEvent;
class KoColorProfile;

/**
 * KisRootSurfaceInfoProxy is a special proxy object for the
 * surface color management information of the widget.
 *
 * When created, the proxy does the following:
 *
 * 1) Finds the top-level widget that \p watched belongs to.
 *    (is \p watched changed its parent in the meantime, e.g.
 *     during construction, then the proxy will handle that
 *     as well)
 *
 * 2) Finds the QWindow that this toplevel widget is painted on
 *    (if platform window is not created yet, subscribes to
 *     QPlatformSurfaceEvent to attach when platform window is
 *     finally created)
 *
 * 3) Finds KisSRGBSurfaceColorSpaceManager object attached to
 *    this window. This object is usually attached with an event
 *    filter in KisApplication to all windows not having a special
 *    property (if the manager is not yet attached, tracks its
 *    addition by filtering QChildEvent)
 *
 * 4) Finally, connects to the manager and forwards its
 *    sigDisplayConfigChanged() to local sigRootSurfaceProfileChanged()
 */
class KRITAUI_EXPORT KisRootSurfaceInfoProxy : public QObject
{
    Q_OBJECT
public:
    KisRootSurfaceInfoProxy(QWidget *watched, QObject *parent = nullptr);
    ~KisRootSurfaceInfoProxy();
    const KoColorProfile* rootSurfaceProfile() const;
    bool isReady() const;

    QString colorManagementReport() const;
    QString osPreferredColorSpaceReport() const;

Q_SIGNALS:
    void sigRootSurfaceProfileChanged(const KoColorProfile *profile) const;

private:

    bool eventFilter(QObject *watched, QEvent *event) override;

    void tryUpdateHierarchy();
    void tryReconnectSurfaceManager();
    void tryUpdateRootSurfaceProfile();

    QVector<QPointer<QObject>> getCurrentHierarchy(QWidget *wdg);

    void reconnectToHierarchy(const QVector<QPointer<QObject>> newHierarchy);

private:
    QWidget *m_watched {nullptr};
    QVector<QPointer<QObject>> m_watchedHierarchy;

    QPointer<QWidget> m_topLevelWidgetWithSurface;

    QPointer<QWindow> m_topLevelNativeWindow;
    QPointer<QObject> m_childChangedFilter;

    QPointer<KisSRGBSurfaceColorSpaceManager> m_topLevelSurfaceManager;
    QMetaObject::Connection m_surfaceManagerConnection;

    const KoColorProfile* m_rootSurfaceProfile {nullptr};

};

#endif /* KISROOTSURFACEINFOPROXY_H */
