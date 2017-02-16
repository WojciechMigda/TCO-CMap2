#!/bin/bash

./run_conc ../data/offline_query_up_n250.csv ../data/offline_query_down_n250.csv ../data/o_wtks_conc.bin #&

exit

PID=$!
chrt -a -b -p 0 $PID
chrt -p $PID
wait $PID
