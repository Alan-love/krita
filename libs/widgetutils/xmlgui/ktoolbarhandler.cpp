/* This file is part of the KDE libraries
   SPDX-FileCopyrightText: 2002 Simon Hausmann <hausmann@kde.org>

   SPDX-License-Identifier: LGPL-2.0-only
*/

#include "ktoolbarhandler_p.h"
#include "config-xmlgui.h"
#include <QDomDocument>
#include <QAction>
#include <QMenu>
#include <QPointer>

#include <kactionmenu.h>
#include <klocalizedstring.h>

#include "kxmlguiwindow.h"
#include "ktoggletoolbaraction.h"
#include "ktoolbar.h"
#include "kxmlguifactory.h"
#include "kactioncollection.h"

#if defined(KCONFIG_BEFORE_5_24)
# define authorizeAction authorizeKAction
#endif

namespace
{
const char actionListName[] = "show_menu_and_toolbar_actionlist";

const char guiDescription[] = ""
                             "<!DOCTYPE kpartgui><kpartgui name=\"StandardToolBarMenuHandler\">"
                             "<MenuBar>"
                             "    <Menu name=\"settings\">"
                             "        <ActionList name=\"%1\" />"
                             "    </Menu>"
                             "</MenuBar>"
                             "</kpartgui>";

class BarActionBuilder
{
public:
    BarActionBuilder(KisKActionCollection *actionCollection, KXmlGuiWindow *mainWindow,
                     QLinkedList<KisToolBar *> &oldToolBarList)
        : m_actionCollection(actionCollection), m_mainWindow(mainWindow), m_needsRebuild(false)
    {
        QList<KisToolBar *> toolBars = m_mainWindow->findChildren<KisToolBar *>();

        Q_FOREACH (KisToolBar *toolBar, toolBars) {
            if (toolBar->mainWindow() != m_mainWindow) {
                continue;
            }

            if (!oldToolBarList.contains(toolBar)) {
                m_needsRebuild = true;
            }

            m_toolBars.append(toolBar);
        }

        if (!m_needsRebuild) {
            m_needsRebuild = (oldToolBarList.count() != m_toolBars.count());
        }
    }

    bool needsRebuild() const
    {
        return m_needsRebuild;
    }

    QList<QAction *> create()
    {
        QList<QAction *> actions;

        if (!m_needsRebuild) {
            return actions;
        }

        Q_FOREACH (KisToolBar *bar, m_toolBars) {
            handleToolBar(bar);
        }

        if (m_toolBarActions.count() == 0) {
            return actions;
        }

        if (m_toolBarActions.count() == 1) {
            KToggleToolBarAction *action = static_cast<KToggleToolBarAction *>(m_toolBarActions.first());
            action->setText(KStandardShortcut::label(KStandardShortcut::ShowToolbar));
            return m_toolBarActions;
        }

        KActionMenu *menuAction = new KActionMenu(i18n("Toolbars Shown"), m_actionCollection);
        m_actionCollection->addAction(QStringLiteral("toolbars_submenu_action"), menuAction);

        Q_FOREACH (QAction *action, m_toolBarActions) {
            menuAction->menu()->addAction(action);
        }

        actions.append(menuAction);

        return actions;
    }

    const QLinkedList<KisToolBar *> &toolBars() const
    {
        return m_toolBars;
    }

private:
    void handleToolBar(KisToolBar *toolBar)
    {
        KToggleToolBarAction *action = new KToggleToolBarAction(
            toolBar,
            toolBar->windowTitle(),
            m_actionCollection);
        m_actionCollection->addAction(toolBar->objectName(), action);

        // ## tooltips, whatsthis?
        m_toolBarActions.append(action);
    }

    KisKActionCollection *m_actionCollection;
    KXmlGuiWindow *m_mainWindow;

    QLinkedList<KisToolBar *> m_toolBars;
    QList<QAction *> m_toolBarActions;

    bool m_needsRebuild : 1;
};
}

using namespace KDEPrivate;

class ToolBarHandler::Private
{
public:
    Private(ToolBarHandler *_parent)
        : parent(_parent)
    {
    }

    void clientAdded(KisKXMLGUIClient *client)
    {
        Q_UNUSED(client);
        parent->setupActions();
    }

    void init(KXmlGuiWindow *mainWindow);
    void connectToActionContainers();
    void connectToActionContainer(QAction *action);
    void connectToActionContainer(QWidget *container);

    ToolBarHandler *parent;
    QPointer<KXmlGuiWindow> mainWindow;
    QList<QAction *> actions;
    QLinkedList<KisToolBar *> toolBars;
};

void ToolBarHandler::Private::init(KXmlGuiWindow *mw)
{
    mainWindow = mw;

    QObject::connect(mainWindow->guiFactory(), SIGNAL(clientAdded(KisKXMLGUIClient*)),
                     parent, SLOT(clientAdded(KisKXMLGUIClient*)));

    if (parent->domDocument().documentElement().isNull()) {

        QString completeDescription = QString::fromLatin1(guiDescription)
                                      .arg(QLatin1String(actionListName));

        parent->setXML(completeDescription, false /*merge*/);
    }
}

void ToolBarHandler::Private::connectToActionContainers()
{
    Q_FOREACH (QAction *action, actions) {
        connectToActionContainer(action);
    }
}

void ToolBarHandler::Private::connectToActionContainer(QAction *action)
{
    uint containerCount = action->associatedWidgets().count();

    for (uint i = 0; i < containerCount; ++i) {
        connectToActionContainer(action->associatedWidgets().value(i));
    }
}

void ToolBarHandler::Private::connectToActionContainer(QWidget *container)
{
    QMenu *popupMenu = qobject_cast<QMenu *>(container);
    if (!popupMenu) {
        return;
    }

    connect(popupMenu, SIGNAL(aboutToShow()),
            parent, SLOT(setupActions()));
}

ToolBarHandler::ToolBarHandler(KXmlGuiWindow *mainWindow)
    : QObject(mainWindow), KisKXMLGUIClient(mainWindow),
      d(new Private(this))
{
    d->init(mainWindow);
}

ToolBarHandler::ToolBarHandler(KXmlGuiWindow *mainWindow, QObject *parent)
    : QObject(parent), KisKXMLGUIClient(mainWindow),
      d(new Private(this))
{
    d->init(mainWindow);
}

ToolBarHandler::~ToolBarHandler()
{
    qDeleteAll(d->actions);
    d->actions.clear();

    delete d;
}

QAction *ToolBarHandler::toolBarMenuAction()
{
    Q_ASSERT(d->actions.count() == 1);
    return d->actions.first();
}

void ToolBarHandler::setupActions()
{
    if (!factory() || !d->mainWindow) {
        return;
    }

    BarActionBuilder builder(actionCollection(), d->mainWindow, d->toolBars);

    if (!builder.needsRebuild()) {
        return;
    }

    unplugActionList(QLatin1String(actionListName));

    qDeleteAll(d->actions);
    d->actions.clear();

    d->actions = builder.create();

    d->toolBars = builder.toolBars();

    // We have no XML file associated with our action collection, so load settings from KConfig
    actionCollection()->readSettings(); // #233712

    plugActionList(QLatin1String(actionListName), d->actions);

    d->connectToActionContainers();
}

#include "moc_ktoolbarhandler_p.cpp"
