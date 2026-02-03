#ifndef SYS_TYPES_H
#define SYS_TYPES_H

// Basic types
typedef unsigned char BYTE;
typedef bool BIT;

#define GET_VALUE_BIT(x, b) (((x) >> (b)) & 1)

/*
struct NormalDistributionParameter {
    float expectation;
    float standard_deviation;
};
*/

enum CellType {
    SLC_CELL = 1,
    MLC_CELL = 2,
    TLC_CELL = 3,
    QLC_CELL = 4,
};

enum OperationType {
    PROGRAM_ONE_SHOT,
    READ_PAGE,
    ERASE_BLOCK,
};

#endif // SYS_TYPES_H
