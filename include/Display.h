/**
 * @file Display.h
 * @brief Handles the curses window thing
*/

/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/
#pragma once

#include <curses.h> // key detection linux
#include <assert.h>
#include <string>
#include <sstream>
//#include <conio.h>   // key detection for windows 

/**
 * @brief   // Curses implementation..
 * @details  This is purely to allow for pausing, and ending, the simualation with a key press. QT would be a more long term solution.
 * @note Keep an eye on effect on computational time. Seems fine currently.
 * https://tldp.org/HOWTO/NCURSES-Programming-HOWTO/keys.html
 */

struct Display{ 
    static void create_screen();
    static void show(const std::stringstream& str);
    static void show(const std::stringstream& str1,const std::stringstream& str2);
    static void clean();
    static void deactivate();
    static void reactivate();
    static void close();
    static constexpr double WIDTH = 30;
    static constexpr double HEIGHT = 10;
    static WINDOW* win;      
    
    static std::string header; // displayed text at start of screen/terminal that doesn't change.

    static std::stringstream display_stream;
    static std::stringstream popup_stream;

    // Could just use this rather than do the screen thing (oops)
    static void signalHandler( int signum );
};


