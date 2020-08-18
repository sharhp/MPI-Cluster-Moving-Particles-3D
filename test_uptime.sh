 for i in 14 11 12 15 16 26 27 28 29 30 42 43 44 46 56 58 59 60 72 73 74 75 76 87 88 89 90 103 104 106 107 117 118 119 120
 do
	echo $i
    ssh csews$i uptime
    ping -c 1 csews$i

 done

