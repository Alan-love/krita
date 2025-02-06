/*
 *  SPDX-FileCopyrightText: 2024 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */
#ifndef KOWRITINGSYSTEMUTILS_H
#define KOWRITINGSYSTEMUTILS_H

#include <QFontDatabase>
#include <QLocale>
#include <kritaflake_export.h>
/**
 * @brief The KoScriptUtils class
 *
 * Collection of utility functions to wrangle the different
 * script and writing system enums in QFontDataBase, QLocale and QChar and ISO 15924 tags
 */

class KRITAFLAKE_EXPORT KoWritingSystemUtils
{
public:
    static QString scriptTagForWritingSystem(QFontDatabase::WritingSystem system);
    static QFontDatabase::WritingSystem writingSystemForScriptTag(const QString &tag);

    // Qt6 has a function to get the ISO 15924 for the QLocale::Script, but we're not qt 6 yet...
    static QString scriptTagForQLocaleScript(QLocale::Script script);
    static QLocale::Script scriptForScriptTag(const QString &tag);

    static QString scriptTagForQCharScript(QChar::Script script);
    static QChar::Script qCharScriptForScriptTag(const QString &tag);

    /**
     * This returns a map of samples and an associated tag. Note that the Sample is the first entry, the tag the second.
     * This is because Latin, for example, has multiple sample strings associated depending on Latin coverage in the font.
     * String it is stored with is s_<ISO 15924> tag for scripts and l_<BCP 47 Language> tag for languages.
     * This way we can have samples per language as is useful for vietnamese.
     */
    static QMap<QString, QString> samples();

    static QString sampleTagForQLocale(const QLocale &locale);
};

#endif // KOWRITINGSYSTEMUTILS_H
