#!/bin/sh
grep ">" $1 | tr -s ' ' | cut -d ' ' -f 1 | sed 's/>//' 1> $2
