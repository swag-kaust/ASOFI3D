#ifndef ENUM_H
#define ENUM_H
/*
 * Enumerations of possible values of integer variables.
 * Allows to avoid magic constants in the code.
 */

// Possible values of the SOURCE_TYPE global variable.
enum SOURCE_TYPE_ENUM {
    // The random type is used for dot-product test.
    SOURCE_TYPE_RANDOM = 0,
    SOURCE_TYPE_EXPLOSIVE = 1,
    SOURCE_TYPE_FORCE_IN_X = 2,
    SOURCE_TYPE_FORCE_IN_Y = 3,
    SOURCE_TYPE_FORCE_IN_Z = 4,
    SOURCE_TYPE_CUSTOM = 5,
    SOURCE_TYPE_EARTHQUAKE = 6,
    SOURCE_TYPE_MOMENT_TENSOR = 7,
};

/*
 * Possible file formats for seismograms and wavefield snapshots.
 *
 * IEEE754/LE format means that the float values are stored according
 * to the IEEE 754 standard [1] for floating-point arithmetic
 * and the bytes are stored in the Little Endian byte order [2].
 * This format is native for commonly used personal computers.
 *
 * IBM/BE format means that the float values are stored according
 * to the IBM standard [3] for floating-point arithmetic
 * and the bytes order is Big Endian [2].
 * This format is "native" for SEG-Y file format.
 *
 * References
 * ----------
 * [1] https://en.wikipedia.org/wiki/IEEE_754
 * [2] https://en.wikipedia.org/wiki/Endianness
 * [3] https://en.wikipedia.org/wiki/IBM_hexadecimal_floating_point
 */
enum FILE_FORMAT_ENUM {
    // SU: Seismic Unix with 4-byte floats in IEEE754/LE format.
    FILE_FORMAT_SU = 1,
    // ASCII: textual format (inefficient for large files, for debug purposes).
    FILE_FORMAT_ASCII = 2,
    // Plain binary file with 4-byte floats in IEEE754/LE format.
    FILE_FORMAT_BINARY = 3,
    // SEG-Y format with 4-byte floats in IEEE754/LE format.
    FILE_FORMAT_SEGY_IEEE754_LITTLEEND = 4,
    // SEG-Y format with 4-byte floats in IBM/BE format.
    FILE_FORMAT_SEGY_IBM_BIGEND = 5,
};
#endif
