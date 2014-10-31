#!/usr/bin/awk -f
BEGIN {
    FS="\t" # use TAB as input field separator
    OFS="\t" # use TAB as output field separator
}
{
    if ($10!=$14 && $9=="+") {
        print
    }
}
