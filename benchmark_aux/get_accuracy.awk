#!/usr/bin/env awk -f
{
    if (NR == 1) {
        printf $1
    }
}
