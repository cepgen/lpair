#!/bin/sh
TEST=`bjobs | grep "$1" | cut -d\  -f1 | awk '/--/{print s;printf $0;s=""}!/--/{s=s" "$0}END{print s}'`
bkill $TEST