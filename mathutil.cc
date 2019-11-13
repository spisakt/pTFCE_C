#include "mathutil.h"

bool isnan2(double a) {
    bool r;
    asm(".intel_syntax noprefix"
        "\n\t ucomisd %1, %1"
        "\n\t setp %b0"
        "\n\t .att_syntax prefix"
        : "=g" (r)
        : "x" (a)
        : "cc"
        );
    return r;
}
