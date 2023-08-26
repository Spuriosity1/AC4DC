/**
 * @file Display.cpp
 * @brief 
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

#include "Display.h"
#include <csignal>
#include <iostream>
#ifdef NCURSES
#include <curses.h> // key detection linux (i.e. TUI).
#endif

// Initialise static variables.
#ifdef NCURSES
WINDOW* Display::win;
#else 
bool Display::win; // Dummy type (bool) that does not use ncurses library.
#endif //NCURSES
std::string Display::header; 
std::stringstream Display::display_stream, Display::popup_stream;

void Display::create_screen(){
    #ifdef NCURSES
    initscr();
    clear();
    noecho();    
    cbreak();
    int startx = (80 - WIDTH) / 2;
    int starty = (24 - HEIGHT) / 2;
    signal(SIGINT,Display::signalHandler);  // clean up for interrupt
    //win = newwin(HEIGHT, WIDTH, starty, startx);
    win = newwin(0, 0, 0, 0);
    box(win, 0 , 0);	
    wrefresh(win);
    keypad(win, TRUE);
    nodelay(win,TRUE); // don't wait for input  
    #endif //NCURSES
}

/// Screen displays the contents of the stream, and only the contents.
void Display::show(const std::stringstream& spooky_stream){
    #ifdef NCURSES
    werase(win); 
    box(win, 0 , 0);
    waddstr(win,spooky_stream.str().c_str());
    wrefresh(win);   
    #endif //NCURSES
}
void Display::show(const std::stringstream& spooky_stream,const std::stringstream& second_stream){
    #ifdef NCURSES
    werase(win);  
    box(win, 0 , 0);
    waddstr(win,(spooky_stream.str()+second_stream.str()).c_str());
    wrefresh(win);   
    #endif //NCURSES
}  
void Display::close(){
    #ifdef NCURSES
    clrtoeol();
    refresh();
    endwin();     
    #endif //NCURSES
}  

// cleans up the terminal on interrupt
void Display::signalHandler( int signum ) {
    #ifdef NCURSES
    endwin();
    std::cout << "Window ended successfully after interrupt signal (" << signum << ") received.\n";
    std::exit(signum);  
    #endif //NCURSES
}    
