//
// Created by Shixuan Sun on 2018/6/29.
//

#ifndef SUBGRAPHMATCHING_COMMANDPARSER_H
#define SUBGRAPHMATCHING_COMMANDPARSER_H

#include <string>
#include <algorithm>
#include <vector>
class CommandParser {
private:
    std::vector<std::string> tokens_;

public:
    CommandParser(const int argc, char **argv);
    const std::string getCommandOption(const std::string &option) const;
    bool commandOptionExists(const std::string &option) const;
};


#endif //SUBGRAPHMATCHING_COMMANDPARSER_H
