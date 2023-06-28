
from IL_3_57_loops_and_strands import loops_and_strands

loopdata = loops_and_strands

for loop_id in loopdata:
	all_nontriv_symmetries = True
	for position, unit_id in loopdata[loop_id]:
		fields = unit_id.split("|")
		if len(fields) < 9 or fields[8] == '1_555':
			all_nontriv_symmetries = False
	if all_nontriv_symmetries:
		del loopdata[loop_id]
		print("Removing %s with %s" % (loop_id,unit_id))

		