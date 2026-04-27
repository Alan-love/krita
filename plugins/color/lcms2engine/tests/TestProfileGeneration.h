/*
 *  SPDX-FileCopyrightText: 2021 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */
#ifndef TESTPROFILEGENERATION_H
#define TESTPROFILEGENERATION_H

#include <QObject>

class TestProfileGeneration : public QObject
{
    Q_OBJECT
private Q_SLOTS:
    void testTransferFunctions();
    void testCICPwriting();

};

#endif // TESTPROFILEGENERATION_H
