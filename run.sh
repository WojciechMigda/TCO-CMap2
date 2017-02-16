#!/bin/bash

./run_conc ${@} #&

exit

PID=$!
chrt -a -b -p 0 $PID
chrt -p $PID
wait $PID
