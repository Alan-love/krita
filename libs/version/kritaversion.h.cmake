/* This file is part of the Krita libraries
    SPDX-FileCopyrightText: 2003 David Faure <faure@kde.org>
    SPDX-FileCopyrightText: 2003 Lukas Tinkl <lukas@kde.org>
    SPDX-FileCopyrightText: 2004 Nicolas Goutte <goutte@kde.org>
    SPDX-FileCopyrightText: 2015 Jarosław Staniek <staniek@kde.org>

    SPDX-License-Identifier: LGPL-2.0-or-later
*/

#ifndef _KRITA_VERSION_H_
#define _KRITA_VERSION_H_

#include "kritaversion_export.h"

// -- WARNING: do not edit values below, instead edit KRITA_* in /CMakeLists.txt --

/**
* @def KRITA_VERSION_STRING
* @ingroup KritaMacros
* @brief Version of Krita as string, at compile time
*
* This macro contains the Krita version in string form. As it is a macro,
* it contains the version at compile time.
*
* @note The version string might contain spaces and special characters,
* especially for development versions of Krita.
* If you use that macro directly for a file format (e.g. OASIS Open Document)
* or for a protocol (e.g. http) be careful that it is appropriate.
* (Fictional) example: "3.0 Alpha"
*/
#define KRITA_VERSION_STRING "@KRITA_VERSION_STRING@"

/**
 * @def KRITA_STABLE_VERSION_MAJOR
 * @ingroup KritaMacros
 * @brief Major version of stable Krita, at compile time
 * KRITA_VERSION_MAJOR is computed based on this value.
*/
#define KRITA_STABLE_VERSION_MAJOR @KRITA_STABLE_VERSION_MAJOR@

/**
 * @def KRITA_STABLE_VERSION_MINOR
 * @ingroup KritaMacros
 * @brief Minor version of stable Krita, at compile time
 * KRITA_VERSION_MINOR is computed based on this value.
 */
#define KRITA_STABLE_VERSION_MINOR @KRITA_STABLE_VERSION_MINOR@

/**
 * @def KRITA_VERSION_RELEASE
 * @ingroup KritaMacros
 * @brief Release version of Krita, at compile time.
 * 89 for Alpha.
 */
#define KRITA_VERSION_RELEASE @KRITA_VERSION_RELEASE@

/**
 * @def KRITA_ALPHA
 * @ingroup KritaMacros
 * @brief If defined (1..9), indicates at compile time that Krita is in alpha stage
 */
#cmakedefine KRITA_ALPHA @KRITA_ALPHA@

/**
 * @def KRITA_BETA
 * @ingroup KritaMacros
 * @brief If defined (1..9), indicates at compile time that Krita is in beta stage
 */
#cmakedefine KRITA_BETA @KRITA_BETA@

/**
 * @def KRITA_RC
 * @ingroup KritaMacros
 * @brief If defined (1..9), indicates at compile time that Krita is in "release candidate" stage
 */
#cmakedefine KRITA_RC @KRITA_RC@

/**
 * @def KRITA_STABLE
 * @ingroup KritaMacros
 * @brief If defined, indicates at compile time that Krita is in stable stage
 */
#cmakedefine KRITA_STABLE @KRITA_STABLE@

/**
 * @ingroup KritaMacros
 * @brief Make a number from the major, minor and release number of a Krita version
 *
 * This function can be used for preprocessing when KRITA_IS_VERSION is not
 * appropriate.
 */
#define KRITA_MAKE_VERSION( a,b,c ) (((a) << 16) | ((b) << 8) | (c))

/**
 * @ingroup KritaMacros
 * @brief Version of Krita as number, at compile time
 *
 * This macro contains the Krita version in number form. As it is a macro,
 * it contains the version at compile time. See version() if you need
 * the Krita version used at runtime.
 */
#define KRITA_VERSION \
    KRITA_MAKE_VERSION(KRITA_STABLE_VERSION_MAJOR,KRITA_STABLE_VERSION_MINOR,KRITA_VERSION_RELEASE)

#endif // _KRITA_VERSION_H_
