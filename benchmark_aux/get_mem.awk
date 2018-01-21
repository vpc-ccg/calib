#!/bin/awk -f
{
if ($0 ~/Maximum resident set size/) {
    printf substr($0, index($0, "Maximum resident set size (kbytes): ") + length("Maximum resident set size (kbytes): "))
}
}
