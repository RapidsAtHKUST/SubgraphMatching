set(UTILITY_SRC
        graphoperations.cpp
        commandparser.cpp
        commandparser.h
        computesetintersection.h
        computesetintersection.cpp
        bitsetoperation.cpp
        bitsetoperation.h
        QFliter.h
        han/intersection_algos.cpp
        han/utils/util.cpp)

add_library(utility SHARED
        ${UTILITY_SRC})