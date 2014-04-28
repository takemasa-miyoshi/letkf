#!/bin/sh
# =======================================================================
#
#	Utility Shell Finctions for WRF_LETKF
#
#                                                   2010.05.11 M.Kunii
# =======================================================================

# -----------------------------
#    cal_date2min
# -----------------------------
cal_date2min () {
(

	# calculate total minutes from 1990.01.01.00UTC

	if [ $# -lt 5 ]; then
		echo "Usage : cal_date2min"
		echo "    cal_date2min [yy] [mm] [dd] [hh] [mn]"
		exit
	fi

	tyy=$1
	tmm=$2
	tdd=$3
	thh=$4
	tmn=$5
	
	total_min=0

	# --- year ---	
	yy=1990
	while [ ${yy} -lt ${tyy} ]
	do
		tmp1=`expr ${yy} % 4` || test 1 -eq 1
		tmp2=`expr ${yy} % 100`	 || test 1 -eq 1
	
		if [ ${tmp1} -eq 0 -a ${tmp2} -ne 0 ]; then
			days=366
		else
			days=365
		fi
		
		total_min=`expr ${total_min} + ${days} \* 1440` || test 1 -eq 1
		yy=`expr ${yy} + 1`	 || test 1 -eq 1
	done

	# --- month ---
	mm=1
	while [ ${mm} -lt ${tmm} ]
	do
		tmp1=`expr ${yy} % 4` || test 1 -eq 1
		tmp2=`expr ${yy} % 100` || test 1 -eq 1
	
		if [ ${tmp1} -eq 0 -a ${tmp2} -ne 0 ]; then
			flag=1
		else
			flag=0
		fi
			
		case ${mm} in
			1|3|5|7|8|10|12)
				days=31
				;;
			4|6|9|11)
				days=30
				;;
			2)
				days=`expr 28 + ${flag}`
				;;
		esac
	
		total_min=`expr ${total_min} + ${days} \* 1440`
		mm=`expr ${mm} + 1`
	done

	# --- day ---
	dd=1
	while [ ${dd} -lt ${tdd} ]
	do
		total_min=`expr ${total_min} + 1440`
		dd=`expr ${dd} + 1`
	done
	
	# --- hour ---	
	hh=0
	while [ ${hh} -lt ${thh} ]
	do		
		total_min=`expr ${total_min} + 60`
		hh=`expr ${hh} + 1`
	done

	# --- minute ---	
	total_min=`expr ${total_min} + ${tmn}`
	
	echo ${total_min}
	
)
}


# -----------------------------
#    cal_min2date
# -----------------------------
cal_min2date () {
(

	# calculate date from total minutes from 1990.01.01.00UTC

	if [ $# -lt 1 ]; then
		echo "Usage : cal_min2date"
		echo "    cal_min2date 89411040"
		exit
	fi

	input_min=$1

	# --- year ---	
	yy=1990
	while :
	do
		tmp1=`expr ${yy} % 4` || test 1 -eq 1
		tmp2=`expr ${yy} % 100` || test 1 -eq 1
	
		if [ ${tmp1} -eq 0 -a ${tmp2} -ne 0 ]; then
			days=366
		else
			days=365
		fi
		
		total_min=`expr ${total_min} + ${days} \* 1440`
		if [ ${total_min} -gt ${input_min} ]; then
			total_min=`expr ${total_min} - ${days} \* 1440` || test 1 -eq 1
			break
		fi
		yy=`expr ${yy} + 1`
	done
	tyy=${yy}

	# --- month ---				
	mm=1
	while :
	do
		tmp1=`expr ${yy} % 4` || test 1 -eq 1
		tmp2=`expr ${yy} % 100` || test 1 -eq 1
	
		if [ ${tmp1} -eq 0 -a ${tmp2} -ne 0 ]; then
			flag=1
		else
			flag=0
		fi
			
		case ${mm} in
			1|3|5|7|8|10|12)
				days=31
				;;
			4|6|9|11)
				days=30
				;;
			2)
				days=`expr 28 + ${flag}`
				;;
		esac			
		
		total_min=`expr ${total_min} + ${days} \* 1440`
		if [ ${total_min} -gt ${input_min} ]; then
			total_min=`expr ${total_min} - ${days} \* 1440` || test 1 -eq 1
			break		
		fi				
		mm=`expr ${mm} + 1`	
	done
	if [ ${mm} -lt 10 ]; then
		tmm="0${mm}"
	else
		tmm=${mm}
	fi
	
	# --- day ---	
	dd=1
	while :
	do
		total_min=`expr ${total_min} + 1440`
		if [ ${total_min} -gt ${input_min} ]; then
			total_min=`expr ${total_min} - 1440` || test 1 -eq 1
			break		
		fi		
		dd=`expr ${dd} + 1`		
	done
	if [ ${dd} -lt 10 ]; then
		tdd="0${dd}"
	else
		tdd=${dd}
	fi	

	# --- hour ---
	hh=0
	while :
	do		
		total_min=`expr ${total_min} + 60`
		if [ ${total_min} -gt ${input_min} ]; then
			total_min=`expr ${total_min} - 60` || test 1 -eq 1
			break		
		fi		
		hh=`expr ${hh} + 1`
	done
	if [ ${hh} -lt 10 ]; then
		thh="0${hh}"
	else
		thh=${hh}
	fi

	# --- minute ---		
	mn=`expr ${input_min} - ${total_min}` || test 1 -eq 1
	if [ ${mn} -lt 10 ]; then
		tmn="0${mn}"
	else
		tmn=${mn}
	fi	
	
	echo "${tyy} ${tmm} ${tdd} ${thh} ${tmn}"
		
)
}

# -----------------------------
#    date_edit
# -----------------------------
date_edit () {
(

	if [ $# -lt 6 ]; then
		echo "Usage : date_edit"
		echo "    date_edit [yyyy] [mm] [dd] [hh] [mn] [dt(min)]"
		echo "    ex) date_edit 201005051200 -180"
		exit
	fi
	
	yy=$1
	mm=$2
	dd=$3
	hh=$4
	mn=$5
	dt=$6

	cal_date2min ${yy} ${mm} ${dd} ${hh} ${mn} ${dt} > tmin.$$
	read tmin < tmin.$$
	
	tmin2=`expr ${tmin} + ${dt}`
	
	cal_min2date ${tmin2}
	rm -f tmin.$$	
    
)
}
