/*
 * SPDX-FileCopyrightText: 2018 boud <boud@valdyas.org>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */
#ifndef DUMMYRESOURCE_H
#define DUMMYRESOURCE_H

#include "KoResource.h"
#include <QDebug>
#include <QRandomGenerator64>
#include <KoMD5Generator.h>
#include <KisResourceTypes.h>

class DummyResource : public KoResource {
public:
    DummyResource(const QString &f, const QString &resourceType = ResourceType::PaintOpPresets)
        : KoResource(f)
        , m_resourceType(resourceType)
    {
        QRandomGenerator64 qrg;
        QByteArray ba(1024, '0');
        for (int i = 0; i < 1024 / 8; i+=8) {
            quint64 v = qrg.generate64();
            ba[i] = v;
        }

        m_something = QString::fromUtf8(ba);

        setMD5Sum(KoMD5Generator::generateHash(ba));

        QImage img(512, 512, QImage::Format_RGB32);
        img.fill(Qt::blue);
        setImage(img);

        // add some random metadata to the resource
        addMetaData("test_metadata", expectedMetaData());

        setValid(true);
    }

    DummyResource(const DummyResource &rhs)
        : KoResource(rhs),
          m_something(rhs.m_something)
    {
    }

    KoResourceSP clone() const override
    {
        return KoResourceSP(new DummyResource(*this));
    }

    bool loadFromDevice(QIODevice *dev, KisResourcesInterfaceSP resourcesInterface) override
    {
        Q_UNUSED(resourcesInterface);

        if (!dev->isOpen()) {
            dev->open(QIODevice::ReadOnly);
        }
        setSomething(QString::fromUtf8(dev->readAll()));
        setValid(true);
        return true;
    }

    bool saveToDevice(QIODevice *dev) const override
    {
        if (!dev->isOpen()) {
            dev->open(QIODevice::WriteOnly);
        }
        dev->write(m_something.toUtf8());
        return true;
    }

    void setSomething(const QString &something)
    {
        m_something = something;
        addMetaData("test_metadata", expectedMetaData());
    }

    QString something() const
    {
        return m_something;
    }

    QString expectedMetaData() const 
    {
        return m_something.left(8);
    }

    QPair<QString, QString> resourceType() const override {
        return QPair<QString, QString>(m_resourceType, "");
    }

private:

    QString m_something;
    QString m_resourceType;
};

#endif // DUMMYRESOURCE_H
