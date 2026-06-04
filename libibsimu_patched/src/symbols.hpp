/*! \file symbols.hpp
 *  \brief Latex symbols chart
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */


const FontLib::Symbolname FontLib::symbols[] = 
{
    { "`a",  "\u00E0" },   /* Latin a with grave, à */
    { "´a",  "\u00E1" },   /* Latin a with acute, á */
    { "^a",  "\u00E2" },   /* Latin a with circumflex, â */
    { "~a",  "\u00E3" },   /* Latin a with tilde, ã */
    { "\"a", "\u00E4" },   /* Latin a with diaeresis, ä */
/*{ "?",   "\u00E5" },*/ /* Latin a with ring above, å */
    { "va",  "\u01CE" },   /* Latin a with caron,  */

    { "`A",  "\u00C0" },   /* Latin capital A with grave, À */
    { "´A",  "\u00C1" },   /* Latin capital A with acute, Á */
    { "^A",  "\u00C2" },   /* Latin capital A with circumflex, Â */
    { "~A",  "\u00C3" },   /* Latin capital A with tilde, Ã */
    { "\"A", "\u00C4" },   /* Latin capital A with diaeresis, Ä */
/*{ "?",   "\u00C5" },*/ /* Latin a with ring above, å */
    
    { "`e",  "\u00E8" },  /* Latin e with grave, è */
    { "´e",  "\u00E9" },  /* Latin e with acute, é */
    { "^e",  "\u00EA" },  /* Latin e with circumflex, ê */
    { "~e",  "\u1EBD" },  /* Latin e with tilde, %Gáº½%@ */
    { "\"e", "\u00EB" },  /* Latin e with diaeresis, ë */
    { "ve",  "\u011B" },  /* Latin e with caron,  */
    
    { "`E",  "\u00C8" },  /* Latin capital E with grave, È */
    { "´E",  "\u00C9" },  /* Latin capital E with acute, É */
    { "^E",  "\u00CA" },  /* Latin capital E with circumflex, Ê */
    { "~E",  "\u1EBC" },  /* Latin capital E with tilde, %Gáº¼%@ */
    { "\"E", "\u00CB" },  /* Latin capital E with diaeresis, Ë */
    { "vE",  "\u011A" },  /* Latin capital e with caron,  */
    
    { "`i",  "\u00EC" },  /* Latin i with grave, ì */
    { "´i",  "\u00ED" },  /* Latin i with acute, í */
    { "^i",  "\u00EE" },  /* Latin i with circumflex, î */
    { "~i",  "\u0129" },  /* Latin i with tilde, %GÄ©%@ */
    { "\"i", "\u00EF" },  /* Latin i with diaeresis, ï */
    
    { "`I",  "\u00CC" },  /* Latin capital I with grave, Ì */
    { "´I",  "\u00CD" },  /* Latin capital I with acute, Í */
    { "^I",  "\u00CE" },  /* Latin capital I with circumflex, Î */
    { "~I",  "\u0128" },  /* Latin capital I with tilde, %GÄ¨%@ */
    { "\"I", "\u00CF" },  /* Latin capital I with diaeresis, Ï */
    
    { "`u",  "\u00F9" },  /* Latin u with grave, ù */
    { "´u",  "\u00FA" },  /* Latin u with acute, ú */
    { "^u",  "\u00FB" },  /* Latin u with circumflex, û */
    { "~u",  "\u0169" },  /* Latin u with tilde, %GÅ©%@ */
    { "\"u", "\u00FC" },  /* Latin u with diaeresis, ü */

    { "`U",  "\u00D9" },  /* Latin capital U with grave, Ù */
    { "´U",  "\u00DA" },  /* Latin capital U with acute, Ú */
    { "^U",  "\u00DB" },  /* Latin capital U with circumflex, Û */
    { "~U",  "\u0168" },  /* Latin capital U with tilde, %GÅ¨%@ */
    { "\"U", "\u00DC" },  /* Latin capital U with diaeresis, Ü */
  
    { "`o",  "\u00F2" }, /* Latin o with grave, ò */
    { "´o",  "\u00F3" }, /* Latin o with acute, ó */
    { "^o",  "\u00F4" }, /* Latin o with circumflex, ô */
    { "~o",  "\u00F5" }, /* Latin o with tilde, õ */
    { "\"o", "\u00F6" }, /* Latin o with diaeresis, ö */

    { "`O",  "\u00D2" }, /* Latin capital O with grave, Ò */
    { "´O",  "\u00D3" }, /* Latin capital O with acute, Ó */
    { "^O",  "\u00D4" }, /* Latin capital O with circumflex, Ô */
    { "~O",  "\u00D5" }, /* Latin capital O with tilde, Õ */
    { "\"O", "\u00D6" }, /* Latin capital O with diaeresis, Ö */

    { "cc",  "\u00E7" }, /* Latin c with cedilla,  */
    { "cC",  "\u00C7" }, /* Latin capital C with cedilla */

    { "´y",  "\u00FD" }, /* Latin y with acute, ý */
    { "\"y", "\u00FF" }, /* Latin y with diaeresis, ÿ */

    { "´Y",  "\u00DD" }, /* Latin capital Y with acute, Ý */
    { "\"Y", "\u0178" }, /* Latin capital Y with diaeresis, %GÅ¸%@ */

    { "~n",  "\u00F1" }, /* Latin n with tilde, ñ */

    { "~N",  "\u00D1" }, /* Latin capital N with tilde, Ñ */

    { "´s",  "\u015B" }, /* Latin s with acute, %GÅ›%@ */
    { "cs",  "\u015F" }, /* Latin s with cedilla,  */
    { "vs",  "\u0161" }, /* Latin s with caron,  */

    { "´S",  "\u015A" }, /* Latin capital S with acute, %GÅš%@ */
    { "cS",  "\u015E" }, /* Latin capital S with cedilla,  */
    { "vS",  "\u0160" }, /* Latin capital S with caron,  */

    { "´z",  "\u017A" }, /* Latin z with acute, %GÅº%@ */
    { "vz",  "\u017E" }, /* Latin z with caron,  */

    { "´Z",  "\u0179" }, /* Latin capital Z with acute, %GÅ¹%@ */
    { "vZ",  "\u017D" }, /* Latin capital Z with caron,  */


    /* Standard symbols */
    { "backslash",  "\\" },
    { "lbrace",     "{" },
    { "rbrace",     "}" },
    { "cent",       "\u00A2" },
    { "pounds",     "\u00A3" },
    { "euro",       "\u20AC" },
    { "S",          "\u00A7" },
    { "copyright",  "\u00A9" },
    { "registered", "\u00AE" },
    { "trademark",  "\u2122" },
    { "P",          "\u00B6" },
    { "ae",         "\u00E6" },
    { "AE",         "\u00C6" },
    { "i",          "\u0131" }, /* dotless i */
    { "OE",         "\u0152" },
    { "oe",         "\u0153" },
  
    { "dag",        "\u2020" },
    { "ddag",       "\u2021" },
    { "bullet",     "\u2022" },
    { "dots",       "\u2026" },
    { "fi",         "\uFB01" },
    { "fl",         "\uFB02" },

    /* Capital Greek letters */
    { "Alpha",      "\u0391" },
    { "Beta",       "\u0392" },
    { "Gamma",      "\u0393" },
    { "Delta",      "\u2206" },
    { "Epsilon",    "\u0395" },
    { "Zeta",       "\u0396" },
    { "Eta",        "\u0397" },
    { "Theta",      "\u0398" },
    { "Iota",       "\u0399" },
    { "Kappa",      "\u039A" },
    { "Lambda",     "\u039B" },
    { "Mu",         "\u039C" },
    { "Nu",         "\u039D" },
    { "Xi",         "\u039E" },
    { "Omicron",    "\u039F" },
    { "Pi",         "\u03A0" },
    { "Rho",        "\u03A1" },
    { "Sigma",      "\u03A3" },
    { "Tau",        "\u03A4" },
    { "Upsilon",    "\u03D2" },
    { "Phi",        "\u03A6" },
    { "Chi",        "\u03A7" },
    { "Psi",        "\u03A8" },
    { "Omega",      "\u2126" },
  
    /* Greek letters */
    { "alpha",      "\u03B1" },
    { "beta",       "\u03B2" },
    { "gamma",      "\u03B3" },
    { "delta",      "\u03B4" },
    { "epsilon",    "\u03B5" },
    { "zeta",       "\u03B6" },
    { "eta",        "\u03B7" },
    { "theta",      "\u03B8" },
    { "vartheta",   "\u03D1" },
    { "iota",       "\u03B9" },
    { "kappa",      "\u03BA" },
    { "lambda",     "\u03BB" },
    { "mu",         "\u00B5" },
    { "nu",         "\u03BD" },
    { "xi",         "\u03BE" },
    { "omicron",    "\u03BF" },
    { "pi",         "\u03C0" },
    { "varpi",      "\u03D6" },
    { "rho",        "\u03C1" },
    { "varsigma",   "\u03C2" },
    { "sigma",      "\u03C3" },
    { "tau",        "\u03C4" },
    { "upsilon",    "\u03C5" },
    { "phi",        "\u03C6" },
    { "varphi",     "\u03D5" },
    { "chi",        "\u03C7" },
    { "psi",        "\u03C8" },
    { "omega",      "\u03C9" },
  
    /* Mathematical symbols */
    { "infty",      "\u221E" },
  
    { "approx",     "\u2248" },
    { "ne",         "\u2260" },
    { "neq",        "\u2260" },
    { "equiv",      "\u2261" },
    { "le",         "\u2264" },
    { "leq",        "\u2264" },
    { "ge",         "\u2265" },
    { "geq",        "\u2265" },
    { "cong",       "\u2245" },
    { "propto",     "\u221D" },

    /* Math, logical */
    { "lnot",       "\u00AC" },
    { "neg",        "\u00AC" },
    { "land",       "\u2227" },
    { "lor",        "\u2228" },
    { "cup",        "\u222A" },
    { "cap",        "\u2229" },
    { "sim",        "\u223C" },
  
    /* Math, delimiters */
    { "langle",     "\u2329" },
    { "rangle",     "\u232A" },
  
    /* Operators */
    { "oplus",      "\u2295" },
    { "otimes",     "\u2297" },
    { "times",      "\u00D7" },
    { "minus",      "\u2212" }, /* Math minus, longer than '-' */
    { "cdot",       "\u22C5" },
    { "pm",         "\u00B1" },
    { "div",        "\u00F7" },
    { "nabla",      "\u2207" },
    { "int",        "\u222B" },
    { "sum",        "\u2211" },
    { "prod",       "\u220F" },
    { "partial",    "\u2202" },
  
    /* Logic, groups */
    { "wp",         "\u2118" },
    { "aleph",      "\u2135" },
    { "Im",         "\u2111" },
    { "Re",         "\u211C" },
    { "forall",     "\u2200" },
    { "ni",         "\u2209" },
    { "exists",     "\u2203" },
    { "in",         "\u2208" },
    { "subset",     "\u2282" },
    { "supset",     "\u2283" },
    { "subseteq",   "\u2286" },
    { "supseteq",   "\u2287" },
    { "nothing",    "\u2205" },
  
    /* Misc */
    { "ast",        "\u2217" },
    { "surd",       "\u221A" },
    { "angle",      "\u2220" },
    { "perp",       "\u22A5" },
    { "therefore",  "\u2234" },
    { "lozenge",    "\u25CA" }, /* unfilled diamond */
  
    /* Arrows */
    { "leftarrow",  "\u2190" },
    { "uparrow",    "\u2191" },
    { "rightarrow", "\u2192" },
    { "downarrow",  "\u2193" },
    { "leftrightarrow", "\u2194" },
    { "Leftarrow",  "\u21D0" },
    { "Uparrow",    "\u21D1" },
    { "Rightarrow", "\u21D2" },
    { "Downarrow",  "\u21D3" },
    { "Leftrightarrow", "\u21D4" },
  
    /* Unsorted -- postscript names */
    { "minute",     "\u2032" },
    { "second",     "\u2033" },
    { "fraction",   "\u2044" },
    { "degree",     "\u00B0" },
    { "florin",     "\u0192" }, /* function f */
    { "suchthat",   "\u220B" },
    { "notsubset",  "\u2284" },
    { "spade",      "\u2660" },
    { "club",       "\u2663" },
    { "heart",      "\u2665" },
    { "diamond",    "\u2666" }, /* filled diamond */

    /* End mark */
    { NULL, NULL }
};



















