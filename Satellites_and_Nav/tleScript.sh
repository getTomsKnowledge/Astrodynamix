#!/bin/bash

# curl the TLE data from Celestrak.com:
curl -o TLE_Input.txt https://celestrak.org/NORAD/elements/gp.php?GROUP=gps-ops&FORMAT=tle

# Get the process ID for this task:
# PID=$!

# Wait for read/write of TLE data:
sleep 5

# End the task:
# kill $PID

# Close gracefully :-)
echo $(date)
echo 'Wrote TLE data to "TLE_Input.txt" from NORAD.  Goodbye!'
exit $?
