//
// Created by Shixuan Sun on 2018/6/29.
//

#include "commandparser.h"

CommandParser::CommandParser(const int argc, char **argv) {
    for (int i = 1; i < argc; ++i)
        tokens_.push_back(std::string(argv[i]));
}

const std::string CommandParser::getCommandOption(const std::string &option) const {

    std::vector<std::string>::const_iterator itr;
    itr = find(tokens_.begin(), tokens_.end(), option);
    if (itr != tokens_.end() && ++itr != tokens_.end()) {
        return *itr;
    }
    return "";
}

bool CommandParser::commandOptionExists(const std::string &option) const {
    return find(tokens_.begin(), tokens_.end(), option) != tokens_.end();
}