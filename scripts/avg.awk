#!/bin/awk
# Calculate average, syntax: avg.awk field-number file
BEGIN { field = ARGV[1]; ARGV[1] = "" }
{ tot += $field; count++; totsq = $field*$field }
END { print tot/count }