#ifndef ENUM_H
#define ENUM_H
/*
 * Enumerations of possible values of integer variables.
 * Allows to avoid magic constants in the code.
 */

// Possible values of the SOURCE_TYPE global variable.
enum SOURCE_TYPE_ENUM {
    SOURCE_TYPE_EXPLOSIVE = 1,
    SOURCE_TYPE_FORCE_IN_X = 2,
    SOURCE_TYPE_FORCE_IN_Y = 3,
    SOURCE_TYPE_FORCE_IN_Z = 4,
    SOURCE_TYPE_CUSTOM = 5,
    SOURCE_TYPE_EARTHQUAKE = 6,
    SOURCE_TYPE_RANDOM = 7,
};
#endif
