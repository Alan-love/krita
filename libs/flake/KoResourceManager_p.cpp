/*
   SPDX-FileCopyrightText: 2006, 2011 Boudewijn Rempt (boud@valdyas.org)
   SPDX-FileCopyrightText: 2007, 2010 Thomas Zander <zander@kde.org>
   SPDX-FileCopyrightText: 2008 Carlos Licea <carlos.licea@kdemail.net>
   SPDX-FileCopyrightText: 2011 Jan Hambrecht <jaham@gmx.net>

   SPDX-License-Identifier: LGPL-2.0-or-later
 */
#include "KoResourceManager_p.h"

#include <QVariant>
#include <FlakeDebug.h>

#include "KoShape.h"
#include "kis_assert.h"
#include "kis_debug.h"

void KoResourceManager::slotResourceInternalsChanged(int key)
{
    KIS_SAFE_ASSERT_RECOVER_RETURN(m_resources.contains(key) || m_abstractResources.contains(key));
    notifyDerivedResourcesChanged(key, m_resources[key]);
    notifyDependenciesAboutTargetChange(key, m_resources[key]);
}

void KoResourceManager::slotAbstractResourceChangedExternal(int key, const QVariant &value)
{
    notifyResourceChanged(key, value);
}

void KoResourceManager::setResource(int key, const QVariant &value)
{
    notifyResourceChangeAttempted(key, value);

    KoDerivedResourceConverterSP converter =
        m_derivedResources.value(key, KoDerivedResourceConverterSP());

    KoAbstractCanvasResourceInterfaceSP abstractResource =
            m_abstractResources.value(key, KoAbstractCanvasResourceInterfaceSP());

    if (abstractResource) {
        const QVariant oldValue = abstractResource->value();
        abstractResource->setValue(value);

        if (m_updateMediators.contains(key)) {
            m_updateMediators[key]->connectResource(value);
        }

        if (oldValue != value) {
            notifyResourceChanged(key, value);
        }
    } else if (converter) {
        const int sourceKey = converter->sourceKey();
        const QVariant oldSourceValue = m_resources.value(sourceKey, QVariant());

        bool valueChanged = false;
        const QVariant newSourceValue = converter->writeToSource(value, oldSourceValue, &valueChanged);

        if (valueChanged) {
            notifyResourceChanged(key, value);
        }

        if (oldSourceValue != newSourceValue) {
            m_resources[sourceKey] = newSourceValue;
            notifyResourceChanged(sourceKey, newSourceValue);
        }
    } else if (m_resources.contains(key)) {
        const QVariant oldValue = m_resources.value(key, QVariant());
        m_resources[key] = value;

        if (m_updateMediators.contains(key)) {
            m_updateMediators[key]->connectResource(value);
        }

        if (oldValue != value) {
            notifyResourceChanged(key, value);
        }
    } else {
        m_resources[key] = value;
        if (m_updateMediators.contains(key)) {
            m_updateMediators[key]->connectResource(value);
        }
        notifyResourceChanged(key, value);
    }
}

void KoResourceManager::notifyResourceChanged(int key, const QVariant &value)
{
    Q_EMIT resourceChanged(key, value);
    notifyDerivedResourcesChanged(key, value);
    notifyDependenciesAboutTargetChange(key, value);
}

void KoResourceManager::notifyDerivedResourcesChanged(int key, const QVariant &value)
{
    QMultiHash<int, KoDerivedResourceConverterSP>::const_iterator it = m_derivedFromSource.constFind(key);
    QMultiHash<int, KoDerivedResourceConverterSP>::const_iterator end = m_derivedFromSource.constEnd();

    while (it != end && it.key() == key) {
        KoDerivedResourceConverterSP converter = it.value();

        if (converter->notifySourceChanged(value)) {
            notifyResourceChanged(converter->key(), converter->readFromSource(value));
        }

        it++;
    }
}

void KoResourceManager::notifyResourceChangeAttempted(int key, const QVariant &value)
{
    Q_EMIT resourceChangeAttempted(key, value);
    notifyDerivedResourcesChangeAttempted(key, value);
}

void KoResourceManager::notifyDerivedResourcesChangeAttempted(int key, const QVariant &value)
{
    QMultiHash<int, KoDerivedResourceConverterSP>::const_iterator it = m_derivedFromSource.constFind(key);
    QMultiHash<int, KoDerivedResourceConverterSP>::const_iterator end = m_derivedFromSource.constEnd();

    while (it != end && it.key() == key) {
        KoDerivedResourceConverterSP converter = it.value();
        notifyResourceChangeAttempted(converter->key(), converter->readFromSource(value));
        it++;
    }
}

void KoResourceManager::notifyDependenciesAboutTargetChange(int targetKey, const QVariant &targetValue)
{
    auto it = m_dependencyFromTarget.find(targetKey);
    while (it != m_dependencyFromTarget.end() && it.key() == targetKey) {
        const int sourceKey = it.value()->sourceKey();

        if (hasResource(sourceKey)) {
            QVariant sourceValue = resource(sourceKey);

            notifyResourceChangeAttempted(sourceKey, sourceValue);
            if (it.value()->shouldUpdateSource(sourceValue, targetValue)) {
                notifyResourceChanged(sourceKey, sourceValue);
            }
        }

        ++it;
    }
}

QVariant KoResourceManager::resource(int key) const
{
    KoAbstractCanvasResourceInterfaceSP abstractResource =
            m_abstractResources.value(key, KoAbstractCanvasResourceInterfaceSP());
    if (abstractResource) {
        return abstractResource->value();
    }

    KoDerivedResourceConverterSP converter =
        m_derivedResources.value(key, KoDerivedResourceConverterSP());

    const int realKey = converter ? converter->sourceKey() : key;
    QVariant value = m_resources.value(realKey, QVariant());

    return converter ? converter->readFromSource(value) : value;
}

void KoResourceManager::setResource(int key, const KoColor &color)
{
    QVariant v;
    v.setValue(color);
    setResource(key, v);
}

void KoResourceManager::setResource(int key, KoShape *shape)
{
    QVariant v;
    v.setValue(shape);
    setResource(key, v);
}

void KoResourceManager::setResource(int key, const KoUnit &unit)
{
    QVariant v;
    v.setValue(unit);
    setResource(key, v);
}

KoColor KoResourceManager::koColorResource(int key) const
{
    if (! m_resources.contains(key)) {
        KoColor empty;
        return empty;
    }
    return resource(key).value<KoColor>();
}

KoShape *KoResourceManager::koShapeResource(int key) const
{
    if (! m_resources.contains(key))
        return 0;

    return resource(key).value<KoShape *>();
}


KoUnit KoResourceManager::unitResource(int key) const
{
    return resource(key).value<KoUnit>();
}

bool KoResourceManager::boolResource(int key) const
{
    if (! m_resources.contains(key))
        return false;
    return m_resources[key].toBool();
}

int KoResourceManager::intResource(int key) const
{
    if (! m_resources.contains(key))
        return 0;
    return m_resources[key].toInt();
}

QString KoResourceManager::stringResource(int key) const
{
    if (! m_resources.contains(key)) {
        QString empty;
        return empty;
    }
    return qvariant_cast<QString>(resource(key));
}

QSizeF KoResourceManager::sizeResource(int key) const
{
    if (! m_resources.contains(key)) {
        QSizeF empty;
        return empty;
    }
    return qvariant_cast<QSizeF>(resource(key));
}

bool KoResourceManager::hasResource(int key) const
{
    if (m_abstractResources.contains(key)) return true;

    KoDerivedResourceConverterSP converter =
        m_derivedResources.value(key, KoDerivedResourceConverterSP());

    const int realKey = converter ? converter->sourceKey() : key;
    return m_resources.contains(realKey);
}

void KoResourceManager::clearResource(int key)
{
    // we cannot remove a derived resource
    if (m_derivedResources.contains(key)) return;

    // we cannot remove an abstract resource either
    if (m_abstractResources.contains(key)) return;

    if (m_resources.contains(key)) {
        m_resources.remove(key);
        notifyResourceChanged(key, QVariant());
    }
}

void KoResourceManager::addDerivedResourceConverter(KoDerivedResourceConverterSP converter)
{
    KIS_SAFE_ASSERT_RECOVER_NOOP(!m_derivedResources.contains(converter->key()));

    if (hasAbstractResource(converter->key()))
        qWarning() << "An abstract resource with the same resource ID exists!";

    m_derivedResources.insert(converter->key(), converter);
    m_derivedFromSource.insert(converter->sourceKey(), converter);
}

bool KoResourceManager::hasDerivedResourceConverter(int key)
{
    return m_derivedResources.contains(key);
}

void KoResourceManager::removeDerivedResourceConverter(int key)
{
    Q_ASSERT(m_derivedResources.contains(key));

    KoDerivedResourceConverterSP converter = m_derivedResources.value(key);
    m_derivedResources.remove(key);
    m_derivedFromSource.remove(converter->sourceKey(), converter);
}

void KoResourceManager::addResourceUpdateMediator(KoResourceUpdateMediatorSP mediator)
{
    KIS_SAFE_ASSERT_RECOVER_NOOP(!m_updateMediators.contains(mediator->key()));
    m_updateMediators.insert(mediator->key(), mediator);
    connect(mediator.data(), SIGNAL(sigResourceChanged(int)), SLOT(slotResourceInternalsChanged(int)));
}

bool KoResourceManager::hasResourceUpdateMediator(int key)
{
    return m_updateMediators.contains(key);
}

void KoResourceManager::removeResourceUpdateMediator(int key)
{
    KIS_SAFE_ASSERT_RECOVER_RETURN(m_updateMediators.contains(key));
    m_updateMediators.remove(key);
}

void KoResourceManager::addActiveCanvasResourceDependency(KoActiveCanvasResourceDependencySP dep)
{
    KIS_SAFE_ASSERT_RECOVER_RETURN(!hasActiveCanvasResourceDependency(dep->sourceKey(), dep->targetKey()));

    m_dependencyFromSource.insert(dep->sourceKey(), dep);
    m_dependencyFromTarget.insert(dep->targetKey(), dep);
}

bool KoResourceManager::hasActiveCanvasResourceDependency(int sourceKey, int targetKey) const
{
    auto it = m_dependencyFromSource.find(sourceKey);

    while (it != m_dependencyFromSource.end() && it.key() == sourceKey) {
        if (it.value()->targetKey() == targetKey) {
            return true;
        }
        ++it;
    }

    return false;
}

void KoResourceManager::removeActiveCanvasResourceDependency(int sourceKey, int targetKey)
{
    KIS_SAFE_ASSERT_RECOVER_RETURN(hasActiveCanvasResourceDependency(sourceKey, targetKey));

    {
        auto it = m_dependencyFromSource.find(sourceKey);
        while (it != m_dependencyFromSource.end() && it.key() == sourceKey) {
            if (it.value()->targetKey() == targetKey) {
                it = m_dependencyFromSource.erase(it);
                break;
            } else {
                ++it;
            }
        }
    }

    {
        auto it = m_dependencyFromTarget.find(targetKey);
        while (it != m_dependencyFromTarget.end() && it.key() == targetKey) {
            if (it.value()->sourceKey() == sourceKey) {
                it = m_dependencyFromTarget.erase(it);
                break;
            } else {
                ++it;
            }
        }
    }
}

bool KoResourceManager::hasAbstractResource(int key)
{
    return m_abstractResources.contains(key);
}

void KoResourceManager::removeAbstractResource(int key)
{
    Q_ASSERT(hasAbstractResource(key));

    KoAbstractCanvasResourceInterfaceSP resourceInterface = m_abstractResources.value(key);
    disconnect(resourceInterface.data(), SIGNAL(sigResourceChangedExternal(int, QVariant)),
               this, SLOT(slotAbstractResourceChangedExternal(int, QVariant)));
    m_abstractResources.remove(key);
}

void KoResourceManager::setAbstractResource(KoAbstractCanvasResourceInterfaceSP resourceInterface)
{
    KIS_SAFE_ASSERT_RECOVER_RETURN(resourceInterface);

    if (hasDerivedResourceConverter(resourceInterface->key()))
        qWarning() << "A derived resource converter with the same resource ID exists!";

    const QVariant oldValue = this->resource(resourceInterface->key());

    KoAbstractCanvasResourceInterfaceSP oldResourceInterface =
        m_abstractResources.value(resourceInterface->key());
    if (oldResourceInterface) {
        disconnect(oldResourceInterface.data(), SIGNAL(sigResourceChangedExternal(int, QVariant)),
                   this, SLOT(slotAbstractResourceChangedExternal(int, QVariant)));
    }

    m_abstractResources[resourceInterface->key()] = resourceInterface;

    connect(resourceInterface.data(), SIGNAL(sigResourceChangedExternal(int, QVariant)),
            this, SLOT(slotAbstractResourceChangedExternal(int, const QVariant&)));

    if (oldValue != resourceInterface->value()) {
        notifyResourceChanged(resourceInterface->key(), resourceInterface->value());
    }
}
