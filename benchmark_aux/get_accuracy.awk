#!/bin/awk -f
{
    if (NR == 1) {
        printf $1
    }
}
