/**
 * @file loadingwidget.hpp
 * @brief 
 */

#ifndef LOADINGWIDGET_CXX_HPP
#define LOADINGWIDGET_CXX_HPP

#include <iostream>
#include <string>
#include <vector>

class LoadingWidget{
public:
    LoadingWidget(){
        LoadingWidget("|/-\\");
    }
    LoadingWidget(const std::string& chars){
        for (auto& c : chars){
            std::string entry {c};
            _versions.push_back(entry);
        }
    }
    LoadingWidget(std::vector<std::string>& v){
        _versions = v;
    }
    std::string next(){
        idx = (idx == _versions.size() -1) ? 0 : idx + 1;
        return _versions[idx];
    }
private:
    std::vector<std::string> _versions;
    int idx=0;
};

std::ostream& operator<<(LoadingWidget lw, std::ostream& os){
    os << lw.next();
    return os;
}


#endif