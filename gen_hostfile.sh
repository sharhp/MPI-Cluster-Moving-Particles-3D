#!/bin/bash

while [ -n "$1" ]
do
    case "$1" in
    -b) begin=$2
        shift
        ;;
    -n) num=$2
        shift
        ;;
    -s) step=$2
        shift
        ;;    
    -ov)ov=0
        ;;
    *) 	echo "Incorrect option '$1'"
        exit 1
        ;;
    esac
    shift
done

if [ -f hostfile ]
then
    rm -f hostfile
fi

# Write self hostname to hostfile
self=$(ifconfig | grep 172.* | awk '{print $2}' | awk -F '.' '{print $4}')
    echo "csews$self" >> hostfile

let count=1
#Removed 57, 105, 106 - Connection refused error
#for i in 11 12 14 15 16 26 27 28 29 30 42 43 44 46 56 58 59 60 72 73 74 75 76 87 88 89 90 103 104 107 117 118 119 120
for i in 16 74 107 119 120 117 27 15 44 58 72 75 118 28 30 42 43 46 56 73 87 88 89 90 103 104 11 12 14 26 29 59 60 76
do
    if [ $count -eq $num ]
    then
        break
    fi
    if [ $i -eq $self ]
    then
        let i+=1
        continue
    fi
    #ping -c 1 -w 500 csews$i &> /dev/null && echo csews$i >> hostfile && let count=count+1
    ssh csews$i uptime &> /dev/null && echo csews$i >> hostfile && let count=count+1
done

if [ $count -ne $num ]
then
    echo "Not enough hosts ready"
    exit 1
fi

