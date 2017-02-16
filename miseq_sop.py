# ################################### INFO ##################################### #
# 16S MiSeq SOP PYTHON SCRIPT
# Author: Amir Shams
# Email: amir.shams84@gmail.com
# Manual Details can be found at: https://www.mothur.org/wiki/MiSeq_SOP
# ################################### IMPORTED LIBRARY ######################### #
import sys
import os
import errno
import datetime
import signal
import logging as log
import time
import subprocess
import traceback
import itertools
import argparse
import multiprocessing
import platform
import collections
import shutil
import pandas
import biom
import numpy
import zipfile
# ################################### GLOBALS ################################## #


fastq_extensions = ['.fastq', '.fq']
CHECK_MARK = "OK"
FAILED_MARK = ":("
DEFAULT_OUTPUTDIR = "/MISEQ_SOP_OUTPUTDIR/"
DEFAULT_TESTDIR = "/MISEQ_SOP_TESTDIR/"
DEFAULT_EXECDIR = "/MISEQ_SOP_EXECDIR/"
DEFAULT_TAXDIR = "/MISEQ_SOP_TAXDIR/"
DEFAULT_INPUTDIR = "/MISEQ_SOP_INPUTDIR/"
DEFAULT_PROCESSORS = str(multiprocessing.cpu_count())
DEFAULT_PREFIX = "MISEQ_SOP"
NO_BETA = True
CURRENT_PATH = "./"
TAXONOMY_EXIST = False
DESIGN_EXIST = True
# ################################### OBJECTS ################################## #


class mothur_process:
	def __init__(self, mothur_input_dictionary):
		for var_name, var_value in mothur_input_dictionary.items():
			setattr(self, var_name, var_value)

	def build_mothur_command(self):
		space = " "
		string = ''
		if hasattr(self, 'nohup_in'):
			string += self.nohup_in + space
		string += self.mothur_exec_path + space + '"#set.dir(output=' + self.outputdir + ');' + space
		string += 'set.logfile(name=mothur.' + self.command + '.logfile);' + space
		string += 'set.current(processors=' + self.processors + ');' + space
		string += self.command + '('
		if hasattr(self, 'parameters') and len(self.parameters) > 0:
			for each_element in self.parameters:
				string += each_element
		string += ');"'
		if hasattr(self, 'nohup_out'):
			string += space + self.nohup_out
		if hasattr(self, 'pid_file'):
			string += ' echo $! > ' + self.pid_file
		report(string)
		print string
		return string

	def execute_mothur_command(self):
		exec_dict = {}
		exec_dict = self.execute([self.build_mothur_command()])
		if exec_dict['exitCode'] != 0:
			print "ERROR occurred!!!!"
			return (False, exec_dict)
		else:
			return(True, exec_dict)

	def execute(self, command, ** kwargs):
		assert isinstance(command, list), "Expected 'command' parameter to be a list containing the process/arguments to execute. Got %s of type %s instead" % (command, type(command))
		assert len(command) > 0, "Received empty list of parameters"
		retval = {
					"exitCode": -1,
					"stderr": u"",
					"stdout": u"",
					"execTime": datetime.timedelta(0),
					"command": None,
					"pid": None
				}
		retval["command"] = command
		log.info("::singleProcessExecuter > At %s, executing \"%s\"" % (datetime.datetime.now(), " ".join(command)))
		cwd = kwargs.get("cwd", os.getcwd())
		#user = kwargs.get("user", os.getuid)
		sheel = kwargs.get("shell", True)
		startDatetime = datetime.datetime.now()
		myPopy = subprocess.Popen(command, cwd=cwd, preexec_fn=os.seteuid(os.getuid()), shell=sheel, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		retval["pid"] = myPopy.pid
		log.debug("::singleProcessExecuter > Command \"%s\" got pid %s" % (" ".join(command), myPopy.pid))
		try:
			retval["stdout"], retval["stderr"] = myPopy.communicate()
			myPopy.wait()
		except OSError, osErr:
			log.debug("::singleProcessExecuter > Got %s %s in myPopy.communicate() when trying get output of command %s. It is probably a bug (more info: http://bugs.python.org/issue1731717)" % (osErr, type(osErr), command[0]))
		except Exception, e:
			log.warn("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s" % (type(e), e, " ".join(command)))
			log.debug("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s. Showing traceback:\n%s" % (type(e), e, " ".join(command), traceback.format_exc()))
			raise
		retval["exitCode"] = myPopy.returncode
		retval["execTime"] = datetime.datetime.now() - startDatetime
		return retval
# ################################### EXECUTIONS ############################### #


def execute_functions(function_name, processors, outputdir, thread_type, func_mode, *args):
	list_of_pid = []
	list_of_stderr = []
	list_of_stdout = []
	if thread_type == 'fork':
		threads = int(processors) + 1
		processors = '1'
	elif thread_type == 'multi':
		threads = 2
	for thread in range(1, threads):
		pid_file = 'run.pid' + str(thread)
		stderr_file = 'stderr' + str(thread)
		stdout_file = 'stdout' + str(thread)
		run_pid = outputdir + pid_file
		stderr = outputdir + stderr_file
		stdout = outputdir + stdout_file
		flag = function_name(processors, outputdir, stderr, stdout, run_pid, *args)
		if flag is False:
			sys.exit(2)
		list_of_pid.append(pid_file)
		list_of_stderr.append(stderr_file)
		list_of_stdout.append(stdout_file)
	flag, stderr = process_monitor(list_of_pid, list_of_stderr, list_of_stdout, outputdir, threads, func_mode)
	if flag is False:
		print "Process monitor failed."
	else:
		pass
	return (True, stderr)


def process_monitor(pid_list, stderr_list, stdout_list, outputdir, threads, mode):

	finished_flag = False
	flag_list = {}
	for each_pid in pid_list:
		flag_list[each_pid] = False
	toolbar_width = threads
	sys.stdout.write("[%s]" % (" " * toolbar_width))
	sys.stdout.flush()
	sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['
	while finished_flag is False:
		for pid_file, stderr_file, stdout_file in itertools.izip(pid_list, stderr_list, stdout_list):
			f = open(outputdir + pid_file)
			the_pid = int(f.read().rstrip())
			if pid_exists(the_pid) is False:
				sys.stdout.write("OK")
				sys.stdout.flush()
				flag_list[pid_file] = True
				flag, stderr = validate_execution(stderr_file, stdout_file, outputdir, mode)
				if flag is False:
					sys.stdout.write(":(")
					report("[:()]")
					print "Error in result of this thread: ", str(the_pid)
					error("Error in result of this thread: " + str(the_pid))
					print "All generated threads killed."
					error("All generated threads killed.")
					kill_pid_list(pid_list, outputdir)
					#print stderr
					print "ABORTING!!!"
					report("ABORTING!!!")
					sys.exit(2)
			if False in flag_list.values():
				finished_flag = False
			else:
				finished_flag = True
		time.sleep(1)
	sys.stdout.write("\n")
	report("[OK]")
	return (True, stderr)


def pid_exists(pid):
	if pid < 0:
		return False
	if pid == 0:
		# According to "man 2 kill" PID 0 refers to every process
		# in the process group of the calling process.
		# On certain systems 0 is a valid PID but we have no way
		# to know that in a portable fashion.
		raise ValueError('invalid PID 0')
	try:
		os.kill(pid, 0)
	except OSError as err:
		if err.errno == errno.ESRCH:
			# ESRCH == No such process
			return False
		elif err.errno == errno.EPERM:
			# EPERM clearly means there's a process to deny access to
			return True
		else:
			# According to "man 2 kill" possible error values are
			# (EINVAL, EPERM, ESRCH)
			raise
	else:
		return True


def validate_execution(stderr, stdout, outputdir, mode):
	if mode == 'usearch':
		f = open(outputdir + stdout, 'rU')
		data = f.read()
		#print data
		if '---Fatal error---' in data or 'Invalid command line' in data:
			return (False, data)
		else:
			return(True, data)
	elif mode == 'mothur':
		f = open(outputdir + stdout, 'rU')
		data = f.read()
		if '[ERROR]:' in data:
			return (False, data)
		else:
			return(True, data)
	elif mode == 'mafft':
		f = open(outputdir + stdout, 'rU')
		data = f.read()
		if 'No such file or directory' in data or 'cannot open file' in data:
			return (False, data)
		else:
			return(True, data)


def kill_pid_list(pid_list, outputdir):
	for pid in pid_list:
		f = open(outputdir + pid)
		the_pid = int(f.read().rstrip())
		try:
			os.kill(the_pid, signal.SIGTERM)
		except OSError:
			pass
		print the_pid, "Killed!"
	return True


def execute(command, ** kwargs):
	assert isinstance(command, list), "Expected 'command' parameter to be a list containing the process/arguments to execute. Got %s of type %s instead" % (command, type(command))
	assert len(command) > 0, "Received empty list of parameters"
	retval = {
			"exitCode": -1,
			"stderr": u"",
			"stdout": u"",
			"execTime": datetime.timedelta(0),
			"command": None,
			"pid": None
		}
	retval["command"] = command
	log.info("::singleProcessExecuter > At %s, executing \"%s\"" % (datetime.datetime.now(), " ".join(command)))
	# print("::singleProcessExecuter > At %s, executing \"%s\"" % (datetime.datetime.now(), " ".join(parameter)))
	cwd = kwargs.get("cwd", os.getcwd())
	sheel = kwargs.get("shell", True)
	startDatetime = datetime.datetime.now()
	myPopy = subprocess.Popen(command, cwd=cwd, preexec_fn=os.seteuid(os.getuid()), shell=sheel, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	retval["pid"] = myPopy.pid
	log.debug("::singleProcessExecuter > Command \"%s\" got pid %s" % (" ".join(command), myPopy.pid))
	try:
		retval["stdout"], retval["stderr"] = myPopy.communicate()
		myPopy.wait()
	except OSError, osErr:
		log.debug("::singleProcessExecuter > Got %s %s in myPopy.communicate() when trying get output of command %s. It is probably a bug (more info: http://bugs.python.org/issue1731717)" % (osErr, type(osErr), command[0]))
	except Exception, e:
		log.warn("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s" % (type(e), e, " ".join(command)))
		log.debug("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s. Showing traceback:\n%s" % (type(e), e, " ".join(command), traceback.format_exc()))
		raise
	retval["exitCode"] = myPopy.returncode
	retval["execTime"] = datetime.datetime.now() - startDatetime
	return retval
# ################################### LOGGING & REPORTING ##################### #


def error(report_string):
	f = open(report_file, "a")
	f.write('###############################  ERROR   ###################################\n')
	f.write(report_string)
	f.write("\n")
	f.write('############################################################################\n')
	f.close()


def report(report_string):
	f = open(report_file, "a")
	f.write(report_string.encode('utf8'))
	f.write("\n")
	f.close()
# ################################### UTILITIES ############################## #


def check_it_and_remove_it(filename, noreport=False):
	try:
		os.remove(filename)
		if noreport is False:
			pass
	except OSError:
		pass


def isFileExist(fname):
	# check the exitence/access/path of a file
	if os.path.isfile(fname) and os.path.exists(fname) and os.access(fname, os.R_OK):
		return True
	else:
		return False


def scandirs(path, container, ext_list, mode=None):
	# scan a spath and grab all files by specified extension
	for root, dirs, names in os.walk(path):
		for currentFile in names:
			path, absname, ext = split_file_name(currentFile)
			if mode is None:
			# Looking for exactly the same extension
				if ext in ext_list:
					container.append(os.path.join(root, currentFile))
			elif mode == 'multiple':
				for each_ext in ext_list:
					if ext in each_ext:
						container.append(os.path.join(root, currentFile))
			elif mode == 'partial':
				# when this extension is part of the actual extension of files
				for each_ext in ext_list:
					if each_ext in ext:
						container.append(os.path.join(root, currentFile))
	if len(container) < 1:
		return False
	return True


def split_file_name(file):
	path = os.path.dirname(file) + '/'
	name = os.path.basename(file)
	if '.' in name:
		ext = '.' + '.'.join(name.split('.')[1:])
		absname = name.split('.')[0]
	else:
		ext = 'no_extension'
		absname = name
	return (path, absname, ext)


def relable_fasta(fasta_file, tag, suffix, new_name):
	if tag is None:
		copy_file(fasta_file, new_name)
	else:
		f = open(fasta_file, 'rU')
		string = ''
		c = 1
		for i in f:
			if i[0] == '>':
				string += '>' + tag + str(c).zfill(6) + suffix + '\n'
				c += 1
			else:
				string += i
		hdl = open(new_name, 'w')
		hdl.write(string)
		hdl.close()
		f.close()
		string = ''
	return True


def copy_file(source, destination):
	shutil.copyfile(source, destination)
	return True


def slugify(target_word):
	text = target_word.replace(' ', r'_').replace('\\', r'_').replace('`', r'_').replace('*', r'_').replace('{', r'_').replace('}', r'_').replace('[', r'_').replace(']', r'_').replace('(', r'_').replace(')', r'_').replace('>', r'_').replace('#', r'_').replace('+', r'_').replace('-', r'_').replace('.', r'_').replace('!', r'_').replace('$', r'_').replace("'", r'_').replace('"', r'_').replace('\/', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_')
	return text


def list_to_string(query_list, delimiter=None):
	# convert list to string by specified delimiter
	if delimiter is None:
		delimiter = '\t'
	return delimiter.join(map(str, query_list))


def write_string_down(new_string, file_name):
	f = open(file_name, 'w')
	f.write(new_string)
	f.close()
	return True


def any_file_to_dict_converter_vertical(file_PATH, index_row=None):
	any_DATAFRAME = get_pandas_DATAFRAME(file_PATH)
	any_vertical_LIST = []
	any_vertical_LIST = any_DATAFRAME.columns.values.tolist()
	any_vertical_DICT = {}
	for each_column in any_vertical_LIST:
		each_column_value_list = any_DATAFRAME[each_column].values.tolist()
		any_vertical_DICT[str(each_column)] = map(str, each_column_value_list)
	any_vertical_LIST = map(str, any_vertical_LIST)
	return (any_vertical_DICT, any_vertical_LIST)


def any_file_to_dict_converter_horizontal(file_PATH, index_col=None):
	any_DATAFRAME = get_pandas_DATAFRAME(file_PATH)
	if index_col is None:
		index_col = 0
	any_horizontal_DICT = {}
	any_horizontal_LIST = []
	for index, row in any_DATAFRAME.iterrows():
		row = map(str, row.tolist())
		any_horizontal_DICT[row[index_col]] = row
		any_horizontal_LIST.append(row[index_col])
	return (any_horizontal_DICT, any_horizontal_LIST)


def get_pandas_DATAFRAME(file_PATH):
	extension = get_extension(file_PATH)
	if extension in ['txt', 'tsv', 'csv']:
		return pandas.read_table(file_PATH, low_memory=False, encoding='utf-8', skip_blank_lines=True)
	elif extension in ['xlsx', 'xls']:
		return pandas.read_excel(file_PATH, sheetname=None)
	else:
		print "Unknow extension"
		error("Unknow extension")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	return False


def get_extension(file_PATH):

	return file_PATH.split('.')[-1].lower()


def zip_file_validation_process(zipped_file_PATH, destination_file_PATH):
	#STEP1: check to see if it is a zip file
	if zipfile.is_zipfile(zipped_file_PATH) is False:
		print "THE INPUT FILE SHOULD BE A ZIPPED FOLDER CONTAINS ALL FASTQ FILES"
		error("THE INPUT FILE SHOULD BE A ZIPPED FOLDER CONTAINS ALL FASTQ FILES")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	#STEP2: CHeck the sanity
	zip_handle = zipfile.ZipFile(zipped_file_PATH, 'r')
	if zip_handle.testzip() is not None:
		print "CRC HEADERS ERROR, BAD ZIPPED FILE"
		error("CRC HEADERS ERROR, BAD ZIPPED FILE")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	for zippy in zip_handle.namelist():
		if zippy.startswith('__MACOSX'):
			continue
		zippy_file_name = os.path.basename(zippy)
		if not zippy_file_name:
			continue
		source = zip_handle.open(zippy)
		target = file(os.path.join(destination_file_PATH, zippy_file_name), "wb")
		with source, target:
			shutil.copyfileobj(source, target)
	return True
# ################################### SPECIFIC ############################### #


def parse_mothur_logfile(logfile, keyword, output_container, mode=None):
	if mode == 'file':
		f = open(logfile, 'rU')
		data = f.read().split('\n')
	else:
		data = logfile
	for line in data.split('\n'):
		if keyword in line:
			output_container.append(line)
	if len(output_container) < 1:
		return False
	return True


def data_prep(input_file_directory, files_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	stderr = ''
	# ###################################################################################################################################### #
	flag, stderr = execute_functions(mothur_make_file, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, input_file_directory)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.files']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], files_PATH)
	return True


def merge_relabel(files_PATH, contigs_file_PATH, groups_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	# ###################################################################################################################################### #
	flag, stderr = execute_functions(mothur_make_contigs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, files_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.trim.contigs.fasta']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], contigs_file_PATH)
		scanned_container_list = []
		extension_list = ['.contigs.groups']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], groups_file_PATH)

	# ###################################################################################################################################### #
	return True


def filter_fasta_file(fasta_file_PATH, contigs_groups_file_PATH, filtered_fasta_file_PATH, filtered_group_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	ambig_limit = '0'
	homop_limit = '8'

	median_length = calculate_fasta_median_length(fasta_file_PATH, 'miseq_sop', mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	# ###################################################################################################################################### #
	flag, stderr = execute_functions(mothur_screen_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, contigs_groups_file_PATH, ambig_limit, homop_limit, median_length)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.good.fasta']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], filtered_fasta_file_PATH)

		scanned_container_list = []
		extension_list = ['.good.txt']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], filtered_group_file_PATH)
	# ###################################################################################################################################### #
	return True


def calculate_fasta_median_length(fasta_file, absname, mothur_exec_path, processors, outputdir):
	summary_file = outputdir + absname + '_length_distributions.txt'
	flag, stderr = execute_functions(mothur_summary_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, fasta_file)
	if flag is False:
		print "Execution of mothur_summary_seqs failed!!!"
	mothur_output_container = []
	extension_list = ['.summary']
	flag = scandirs(outputdir, mothur_output_container, extension_list)
	if flag is False:
		print "This extension is not availble: ", extension_list
		sys.exit(2)
	else:
		for each_file in mothur_output_container:
			if 'current_files.summary' in each_file:
				continue
			os.rename(each_file, summary_file)
	the_list = []
	f = open(summary_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == 'seqname':
			continue
		the_list.append(int(line[3]))
	f.close()
	freq_dict = {}
	freq_dict = collections.Counter(the_list)
	freq_list = []
	freq_list = freq_dict.most_common(10)
	top_length, top_freq = grab_most_frequent_length(freq_list)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	freq_string = 'LENGTH FREQUENCY REPORT:\n'
	for i in range(len(freq_list)):
		c = i + 1
		freq_string += 'RANK_' + str(c) + ': ' + str(freq_list[i][0]) + '= ' + str(freq_list[i][1]) + '\n'
	print freq_string
	report(freq_string)
	print "MEDIAN LENGTH:", top_length, " With the frequency of ", top_freq
	report("MEDIAN LENGTH:" + top_length + " With the frequency of " + top_freq)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	return top_length


def grab_most_frequent_length(frequency_list):
	most_frequent_set = frequency_list[0]
	top_length = most_frequent_set[0]
	top_freq = most_frequent_set[1]
	top_limit = int(top_freq / 3)
	for i in range(len(frequency_list)):
		if top_length == frequency_list[i][0]:
			continue
		elif frequency_list[i][0] > top_length:
			continue
		elif frequency_list[i][0] < top_length and top_limit < frequency_list[i][1]:
			top_length = frequency_list[i][0]
			top_freq = frequency_list[i][1]
			top_limit = int(top_freq / 3)
	return (str(top_length), str(top_freq))


def removing_duplicates(fasta_file_PATH, unique_fasta_file_PATH, unique_name_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_unique_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.unique.fasta']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], unique_fasta_file_PATH)

		scanned_container_list = []
		extension_list = ['.names']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], unique_name_file_PATH)
	return True


def counting_sequences(contigs_groups_PATH, unique_names_file_PATH, count_seqs_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_count_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, contigs_groups_PATH, unique_names_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.count_table']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], count_seqs_file_PATH)
	return True


def precluster(fasta_file_PATH, count_file_PATH, precluster_fasta_file_PATH, precluster_count_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_precluster, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, count_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.precluster.fasta']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], precluster_fasta_file_PATH)

		scanned_container_list = []
		extension_list = ['.precluster.count_table']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], precluster_count_file_PATH)
	return True


def chimera_scanning(fasta_file_PATH, count_file_PATH, chimera_report_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_chimera_vsearch, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, count_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.denovo.vsearch.accnos']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], chimera_report_file_PATH)
	return True


def remove_chimeric_seqs(fasta_file_PATH, count_file_PATH, chimera_report_file_PATH, chimera_removed_fasta_file_PATH, chimera_removed_count_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_remove_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, count_file_PATH, chimera_report_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.pick.fasta']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], chimera_removed_fasta_file_PATH)

		scanned_container_list = []
		extension_list = ['.pick.txt']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], chimera_removed_count_file_PATH)
	return True


def classify_seqs(fasta_file_PATH, count_file_PATH, reference_file_PATH, taxonomy_file_PATH, classified_reads_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_classify_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, count_file_PATH, reference_file_PATH, taxonomy_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.knn.taxonomy']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], classified_reads_file_PATH)
	return True


def distance_matrix(fasta_file_PATH, distance_matrix_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	cutoff_limit = '0.03'
	flag, stderr = execute_functions(mothur_dist_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, cutoff_limit)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.dist']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], distance_matrix_file_PATH)
	return True


def cluster_split(distance_matrix_file_PATH, count_table_file_PATH, cluster_list_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_cluster_split, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, distance_matrix_file_PATH, count_table_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.an.unique_list.list']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], cluster_list_file_PATH)
	return True


def filter_list_file(list_file_PATH, filtered_list_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_remove_rare, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, list_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.pick.txt']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], filtered_list_file_PATH)
	return True


def get_representative_OTU(list_file_PATH, count_file_PATH, distance_matrix_file_PATH, fasta_file_PATH, representative_OTU_fasta_file_PATH, representative_OTU_count_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_get_oturep, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, list_file_PATH, count_file_PATH, distance_matrix_file_PATH, fasta_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.rep.fasta']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], representative_OTU_fasta_file_PATH)

		scanned_container_list = []
		extension_list = ['.rep.count_table']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], representative_OTU_count_file_PATH)
	return True


def making_shared(list_file_PATH, count_file_PATH, shared_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_make_shared, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, list_file_PATH, count_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.shared']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], shared_file_PATH)
	return True


def classify_centroids(fasta_file_PATH, reference_file_PATH, taxonomy_file_PATH, classification_report_file_PATH, centroids_classification_fasta_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	
	processed_taxonomy_file_PATH = outputdir_PATH + 'processed_taxonomy_file.txt'
	flag = process_taxonomy_file(taxonomy_file_PATH, processed_taxonomy_file_PATH)
	if flag is True:
		pass

	flag, stderr = execute_functions(mothur_classify_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, reference_file_PATH, processed_taxonomy_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.knn.taxonomy']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		#os.rename(scanned_container_list[0], classification_report_file_PATH)
		uniq_otus_string = ''
		uniq_otu_accnos_file = ''
		f = open(scanned_container_list[0], 'rU')
		uniq_list = []
		for i in f:
			i = i.rstrip()
			line = i.split('\t')
			clean_tax = line[1].split('_TAG_')[0] + ';'
			if clean_tax not in uniq_list:
				uniq_list.append(clean_tax)
				uniq_otus_string += i + '\n'
				uniq_otu_accnos_file += line[0] + '\n'
		f.close()
		write_string_down(uniq_otus_string, classification_report_file_PATH)
		accnos_file = outputdir_PATH + 'accnos_file.txt'
		write_string_down(uniq_otu_accnos_file, accnos_file)
	flag, stderr = execute_functions(mothur_get_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, accnos_file)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.pick.fasta']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], centroids_classification_fasta_file_PATH)
	return True


def global_alignment(fasta_file_PATH, reference_file_PATH, global_alignment_report_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	flag, stderr = execute_functions(mothur_align_seqs, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, fasta_file_PATH, reference_file_PATH)
	if flag is False:
		error("[" + FAILED_MARK + "]")
		print "[", FAILED_MARK, "]"
		print "Execution of mothur failed!!!"
		error("Execution of mothur failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		#it has the file name in it
		scanned_container_list = []
		extension_list = ['.align.report']
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
		os.rename(scanned_container_list[0], global_alignment_report_file_PATH)
	return True


def alignment_report_parser(alignment_report_file_PATH, centroids_taxonomy_file_PATH, groups_file_PATH, EXCEL_file_PATH, SAS_file_PATH, BIOM_file_PATH, shared_file_PATH, outputdir_PATH):
	#STEP1: grab the sample names
	Read_DICT = {}
	#Read_DICT['QueryName'] = 'NO_SAMPLE'
	Sample_name_LIST = []
	f = open(groups_file_PATH, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		Read_DICT[line[0]] = line[1]
	f.close()
	Sample_name_LIST = list(set(Read_DICT.values()))
	#STEP2: FILL UP OTU DICT
	OTU_DICT = {}
	OTU_name_LIST = []
	#OTU_DICT['TemplateName'] = ('ACC', 'NO_TAX;NO_TAX;NO_TAX;NO_TAX;NO_TAX;NO_TAX;NO_TAX;')
	f = open(centroids_taxonomy_file_PATH, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		ACC = line[1].split('_TAG_')[1]
		ACC = ACC[:-1]
		OTU_DICT[line[0]] = (ACC, line[1])
		OTU_name_LIST.append(line[0])
	f.close()

	sas_string = ''
	sas_string += 'Library_id' + '\t' + 'Read_id' + '\t' + 'Query_length' + '\t' + 'Reference_length' + '\t' + 'Compressed_alignment' + '\t' + 'Expected_error_rate' + '\t' + 'Average_error_probability' + '\t' + 'Average_Phred_score' + '\t'
	sas_string += 'Identity' + '\t' + 'Ref_id' + '\t' + 'Ref_acc' + '\t' + 'BEI' + '\t'
	sas_string += 'Kingdom_name' + '\t' + 'Phylum_name' + '\t' + 'Class_name' + '\t' + 'Order_name' + '\t'
	sas_string += 'Family_name' + '\t' + 'Genus_name' + '\t' + 'Species_name' + '\n'

	#OTU_DICT
	satisfy = True
	Sample_DICT = {}
	f = open(alignment_report_file_PATH, 'rU')
	f.next()
	for i in f:
		i = i.rstrip()
		line = i.split('\t')

		QueryName = line[0]
		QueryLength = line[1]
		TemplateName = line[2]
		TemplateLength = line[3]
		SearchMethod = line[4]
		SearchScore = line[5]
		AlignmentMethod = line[6]
		QueryStart = line[7]
		QueryEnd = line[8]
		TemplateStart = line[9]
		TemplateEnd = line[10]
		PairwiseAlignmentLength = line[11]
		GapsInQuery = line[12]
		GapsInTemplate = line[13]
		LongestInsert = line[14]
		SimBtwnQuery_Template = line[15]
		if float(SimBtwnQuery_Template) > 80.0:
			satisfy = True
		else:
			satisfy = False

		######################FLIING SAS STRING
		Library_id = Read_DICT[QueryName]
		Read_id = QueryName
		Query_length = QueryLength
		Reference_length = TemplateLength
		Compressed_alignment = ''
		Expected_error_rate = '.'
		Average_error_probability = '.'
		Average_Phred_score = '.'
		Identity = SimBtwnQuery_Template
		Ref_id = OTU_DICT[TemplateName][0]
		Ref_acc = OTU_DICT[TemplateName][0]
		BEI = '0'
		TAX_LIST = OTU_DICT[TemplateName][1].split(';')
		Kingdom_name = TAX_LIST[0]
		Phylum_name = TAX_LIST[1]
		Class_name = TAX_LIST[2]
		Order_name = TAX_LIST[3]
		Family_name = TAX_LIST[4]
		Genus_name = TAX_LIST[5]
		Species_name = TAX_LIST[6].split('_TAG_')[0]
		if satisfy is True:
			sas_string += Library_id + '\t' + Read_id + '\t' + Query_length + '\t' + Reference_length + '\t' + Compressed_alignment + '\t' + Expected_error_rate + '\t'
			sas_string += Average_error_probability + '\t' + Average_Phred_score + '\t' + Identity + '\t' + Ref_id + '\t' + Ref_acc + '\t' + BEI + '\t'
			sas_string += Kingdom_name + '\t' + Phylum_name + '\t' + Class_name + '\t' + Order_name + '\t' + Family_name + '\t' + Genus_name + '\t' + Species_name + '\n'
		else:
			sas_string += Library_id + '\t' + Read_id + '\t' + Query_length + '\t' + '' + '\t' + '' + '\t' + '.' + '\t'
			sas_string += '.' + '\t' + '.' + '\t' + '.' + '\t' + '' + '\t' + '' + '\t' + '.' + '\t'
			sas_string += '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n'
		######################END OF FILLING SAS STRING
		if satisfy is True:
			Sample_name = Read_DICT[QueryName]
			if Sample_name not in Sample_DICT:
				Sample_DICT[Sample_name] = {}
			elif TemplateName not in Sample_DICT[Sample_name]:
				Sample_DICT[Sample_name][TemplateName] = 1
			else:
				Sample_DICT[Sample_name][TemplateName] += 1

	f.close()

	write_string_down(sas_string, SAS_file_PATH)

	# ##################################SHARED FILE
	sas_string = ''
	numOtus_value = str(len(OTU_name_LIST))
	processed_OTU_name_LIST = []
	for each_OTU in OTU_name_LIST:
		temp_string = 'GG_ID_' + OTU_DICT[each_OTU][0] + ';' + OTU_DICT[each_OTU][1].split('_TAG_')[0] + ';'
		processed_OTU_name_LIST.append(temp_string)
	shared_file_string = 'label\tGroups\tnumOtus\t' + list_to_string(processed_OTU_name_LIST, '\t') + '\n'
	#Sample_name_LIST.remove('NO_SAMPLE')
	for each_Sample in Sample_name_LIST:
		shared_file_string += 'Miseq_SOP' + '\t' + each_Sample + '\t' + numOtus_value + '\t'
		for each_OTU in OTU_name_LIST:
			if each_OTU not in Sample_DICT[each_Sample]:
				shared_file_string += '0' + '\t'
			else:
				shared_file_string += str(Sample_DICT[each_Sample][each_OTU]) + '\t'
		shared_file_string = shared_file_string[:-1]
		shared_file_string += '\n'
	write_string_down(shared_file_string, shared_file_PATH)
	# ##################################END OF SHARED FILE
	# ################################## BIOM FILE GENERATING
	flag = shared_to_biom_converter(shared_file_PATH, OTU_DICT, processed_OTU_name_LIST, BIOM_file_PATH)
	if flag is True:
		pass
	# ################################## END OF BIOM FILE GENERATING
	"""
	temp_Excel_path = outputdir_PATH + 'Raw_Excel_file.xlsx'
	flag = biom_to_excel_converter(BIOM_file_PATH, temp_Excel_path, None)
	if flag is True:
		pass
	flag = interpret_excel_file(temp_Excel_path, EXCEL_file_PATH)
	if flag is True:
		pass
	"""
	return True


def shared_to_biom_converter(shared_file_PATH, OTU_DICT, OTU_name_LIST, BIOM_file_PATH):
	shared_DICT, shared_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	OTU_ID = []
	for each_OTU in OTU_name_LIST:
		temp_ID = each_OTU.split(';')[0]
		OTU_ID.append(temp_ID.split('GG_ID_')[1])

	obs_ids = OTU_ID
	samp_ids = shared_DICT['Groups']
	data_LIST = []
	observe_LIST = []

	for each_OTU in OTU_name_LIST:
		observe_LIST.append({'taxonomy': each_OTU.split(';')[1:-1]})
		data_LIST.append(map(float, shared_DICT[each_OTU]))
	biom_table = biom.table.Table(numpy.array(data_LIST), obs_ids, samp_ids, observe_LIST)
	biom_handle = open(BIOM_file_PATH, 'w')
	biom_table.to_json('miseq_sop', direct_io=biom_handle)
	return True


def biom_to_excel_converter(biom_file_PATH, EXCEL_file_PATH, Design_file_PATH):
	EXCEL_writer = pandas.ExcelWriter(EXCEL_file_PATH)
	#Biom Parsing
	biom_table_HANDLE = biom.load_table(biom_file_PATH)
	
	#DEFINE OUTPUTS
	SHARED_file_DataFrame_DICT = {}
	SHARED_file_DataFrame_header_LIST = []
	#Step1: FILLOUT shared_file_DataFrame_DICT
	SHARED_file_DataFrame_header_LIST = ['label', 'Groups', 'numOtus']
	SHARED_file_DataFrame_DICT['label'] = []
	SHARED_file_DataFrame_DICT['Groups'] = []
	SHARED_file_DataFrame_DICT['numOtus'] = []
	
	sample_ID_LIST = map(str, biom_table_HANDLE.ids(axis='sample'))
	OTU_ID_LIST = map(str, biom_table_HANDLE.ids(axis='observation'))
	
	numOtus_value = len(OTU_ID_LIST)
	alias_counter = 1
	for each_sample_ID in sample_ID_LIST:
		sample_alias_name = 'SAMPLE_ALIAS_' + str(alias_counter).zfill(3)
		SHARED_file_DataFrame_DICT['Groups'].append(sample_alias_name)
		SHARED_file_DataFrame_DICT['label'].append('BIOM_SLICER')
		SHARED_file_DataFrame_DICT['numOtus'].append(numOtus_value)
		alias_counter += 1

	alias_counter = 1
	for each_OTU_ID in OTU_ID_LIST:
		OTU_alias_name = 'OTU_ALIAS_' + str(alias_counter).zfill(6)
		each_OTU_ID_data_LIST = biom_table_HANDLE.data(each_OTU_ID, axis='observation').tolist()
		SHARED_file_DataFrame_DICT[OTU_alias_name] = map(int, each_OTU_ID_data_LIST)
		SHARED_file_DataFrame_header_LIST.append(OTU_alias_name)
		alias_counter += 1
	
	SHARED_DATAFrame = pandas.DataFrame.from_dict(SHARED_file_DataFrame_DICT)
	SHARED_DATAFrame.to_excel(EXCEL_writer, sheet_name='Biome', columns=SHARED_file_DataFrame_header_LIST, index=False, engine='xlsxwriter')
	##################################################################################################################
	##################################################################################################################
	#Step2: FILL OUT SAMPLE_metadata_DataFrame_DICT
	SAMPLE_metadata_DATAFrame_DICT = {}
	SAMPLE_metadata_DATAFrame_header_LIST = []
	SAMPLE_metadata_DATAFrame_DICT['SAMPLE_ALIAS'] = []
	SAMPLE_metadata_DATAFrame_DICT['SAMPLE_ID'] = []
	SAMPLE_metadata_DATAFrame_header_LIST = ['SAMPLE_ALIAS', 'SAMPLE_ID']
	
	sample_ID_LIST = map(str, biom_table_HANDLE.ids(axis='sample'))

	if biom_table_HANDLE.metadata(axis='sample') is not None:
		metadata_sample_OBJECT = biom_table_HANDLE.metadata(axis='sample')
		for each_extra_headers in metadata_sample_OBJECT[0].keys():
			SAMPLE_metadata_DATAFrame_DICT[each_extra_headers] = []
			SAMPLE_metadata_DATAFrame_header_LIST.append(each_extra_headers)

	alias_counter = 1
	for each_sample_ID in sample_ID_LIST:
		sample_alias_name = 'SAMPLE_ALIAS_' + str(alias_counter).zfill(3)
		SAMPLE_metadata_DATAFrame_DICT['SAMPLE_ALIAS'].append(sample_alias_name)
		SAMPLE_metadata_DATAFrame_DICT['SAMPLE_ID'].append(each_sample_ID)

		if biom_table_HANDLE.metadata(each_sample_ID, axis='sample') is not None:
			each_sample_ID_metadata_DICT = dict(biom_table_HANDLE.metadata(each_sample_ID, axis='sample'))
			metadata_keys = each_sample_ID_metadata_DICT.keys()
			for each_metadata_key in metadata_keys:
				if type(each_sample_ID_metadata_DICT[each_metadata_key]) is list:
					SAMPLE_metadata_DATAFrame_DICT[each_metadata_key].append(list_to_string(each_sample_ID_metadata_DICT[each_metadata_key], ';'))
				else:
					SAMPLE_metadata_DATAFrame_DICT[each_metadata_key].append(each_sample_ID_metadata_DICT[each_metadata_key])
		alias_counter += 1
	SAMPLE_metadata_data_frame = pandas.DataFrame.from_dict(SAMPLE_metadata_DATAFrame_DICT)
	SAMPLE_metadata_data_frame.to_excel(EXCEL_writer, sheet_name='Sample_metadata', columns=SAMPLE_metadata_DATAFrame_header_LIST, index=False, engine='xlsxwriter')
	##################################################################################################################
	##################################################################################################################
	#Step#: ENTER THE DESIGN FILE
	if Design_file_PATH is not None:
		DESIGN_metadata_vertical_DICT, DESIGN_metadata_vertical_LIST = any_file_to_dict_converter_vertical(Design_file_PATH)

		key_column = ''
		key_column_len = 0
		#print SAMPLE_metadata_DATAFrame_DICT['SAMPLE_ID']
		for each_Column in DESIGN_metadata_vertical_LIST:
			match_rank = len(match_two_list(SAMPLE_metadata_DATAFrame_DICT['SAMPLE_ID'], DESIGN_metadata_vertical_DICT[each_Column]))
			if match_rank > key_column_len:
				key_column_len = match_rank
				key_column = each_Column
		if key_column_len > 1:
			print "Key Column for design file parsing is: ", key_column
			Category_DICT = {}
			Category_LIST = ['SAMPLE_ALIAS', 'SAMPLE_ID']
			DESIGN_SAMPLE_ALIAS_DICT = {}
			Category_DICT['SAMPLE_ALIAS'] = []
			Category_DICT['SAMPLE_ID'] = []
			for each_SAMPLE_ALIAS, each_SAMPLE_ID in itertools.izip(SAMPLE_metadata_DATAFrame_DICT['SAMPLE_ALIAS'], SAMPLE_metadata_DATAFrame_DICT['SAMPLE_ID']):
				if each_SAMPLE_ID in DESIGN_metadata_vertical_DICT[key_column]:
					DESIGN_SAMPLE_ALIAS_DICT[each_SAMPLE_ID] = each_SAMPLE_ALIAS
			for each_SAMPLE in DESIGN_metadata_vertical_DICT[key_column]:
				Category_DICT['SAMPLE_ALIAS'].append(DESIGN_SAMPLE_ALIAS_DICT[each_SAMPLE])
				Category_DICT['SAMPLE_ID'].append(each_SAMPLE)
			DESIGN_metadata_vertical_LIST.remove(key_column)
			for each_Column in DESIGN_metadata_vertical_LIST:
				Category_LIST.append(each_Column)
				Category_DICT[each_Column] = DESIGN_metadata_vertical_DICT[each_Column]
			DESIGN_metadata_data_frame = pandas.DataFrame.from_dict(Category_DICT)
			DESIGN_metadata_data_frame.to_excel(EXCEL_writer, sheet_name='Design', columns=Category_LIST, index=False, engine='xlsxwriter')
		else:
			print "Design file dose not match with biom file, check the sample name column, will skip parsing design file."
	##################################################################################################################
	##################################################################################################################
	#Step3: FILL OUT OTU_metadata_DataFrame_DICT
	OTU_metadata_DataFrame_DICT = {}
	OTU_metadata_DATAFrame_header_LIST = []
	OTU_metadata_DataFrame_DICT['OTU_ALIAS'] = []
	OTU_metadata_DataFrame_DICT['OTU_ID'] = []
	OTU_metadata_DATAFrame_header_LIST = ['OTU_ALIAS', 'OTU_ID']

	OTU_ID_LIST = map(str, biom_table_HANDLE.ids(axis='observation'))

	if biom_table_HANDLE.metadata(axis='observation') is not None:
		metadata_OTU_OBJECT = biom_table_HANDLE.metadata(axis='observation')
		for each_extra_headers in metadata_OTU_OBJECT[0].keys():
			OTU_metadata_DataFrame_DICT[each_extra_headers] = []
			OTU_metadata_DATAFrame_header_LIST.append(each_extra_headers)

	alias_counter = 1
	for each_OTU_ID in OTU_ID_LIST:
		OTU_alias_name = 'OTU_ALIAS_' + str(alias_counter).zfill(6)
		
		OTU_metadata_DataFrame_DICT['OTU_ALIAS'].append(OTU_alias_name)
		OTU_metadata_DataFrame_DICT['OTU_ID'].append(each_OTU_ID)

		if biom_table_HANDLE.metadata(each_OTU_ID, axis='observation') is not None:
			each_OTU_ID_metadata_DICT = dict(biom_table_HANDLE.metadata(each_OTU_ID, axis='observation'))
			metadata_keys = each_OTU_ID_metadata_DICT.keys()
			for each_metadata_key in metadata_keys:
				if each_metadata_key == 'taxonomy':
					corrected_taxonomy_string = remove_ambiguity_from_taxonomy(each_OTU_ID_metadata_DICT[each_metadata_key])
					OTU_metadata_DataFrame_DICT[each_metadata_key].append(corrected_taxonomy_string)
				elif type(each_OTU_ID_metadata_DICT[each_metadata_key]) is list:
					OTU_metadata_DataFrame_DICT[each_metadata_key].append(list_to_string(each_OTU_ID_metadata_DICT[each_metadata_key], ';'))
				else:
					OTU_metadata_DataFrame_DICT[each_metadata_key].append(each_OTU_ID_metadata_DICT[each_metadata_key])
		alias_counter += 1
	OTU_metadata_data_frame = pandas.DataFrame.from_dict(OTU_metadata_DataFrame_DICT)
	OTU_metadata_data_frame.to_excel(EXCEL_writer, sheet_name='OTU_metadata', columns=OTU_metadata_DATAFrame_header_LIST, index=False, engine='xlsxwriter')
	##################################################################################################################
	##################################################################################################################

	#Step4: FILL OUT MULTI_SHARED SHARED FILE
	lineage_Name_LIST = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
	lineage_target_DICT = {}
	lineage_target_DICT['Kingdom'] = []
	lineage_target_DICT['Phylum'] = []
	lineage_target_DICT['Class'] = []
	lineage_target_DICT['Order'] = []
	lineage_target_DICT['Family'] = []
	lineage_target_DICT['Genus'] = []
	lineage_target_DICT['Species'] = []
	#STEP1: Reduce the number of similar taxonomy
	unique_taxonomy_LIST = []
	unique_taxonomy_LIST = list(set(OTU_metadata_DataFrame_DICT['taxonomy']))
	for each_taxonomy in unique_taxonomy_LIST:
		taxonomy_LIST = each_taxonomy.split(';')
		lineage_target_DICT['Kingdom'].append(taxonomy_LIST[0])
		lineage_target_DICT['Phylum'].append(taxonomy_LIST[1])
		lineage_target_DICT['Class'].append(taxonomy_LIST[2])
		lineage_target_DICT['Order'].append(taxonomy_LIST[3])
		lineage_target_DICT['Family'].append(taxonomy_LIST[4])
		lineage_target_DICT['Genus'].append(taxonomy_LIST[5])
		lineage_target_DICT['Species'].append(taxonomy_LIST[6])
	for each_lineage in lineage_Name_LIST:
		lineage_target_DICT[each_lineage] = list(set(lineage_target_DICT[each_lineage]))
	#print lineage_target_DICT
	#{'Kingdom': ['k_Unclassified', 'k_Archaea', 'k_Bacteria'], 'Family': ['*o_Exiguobacterales', 'f_Fimbriimonadaceae'], 'Species': ['*o_RF39', '*f_Peptococcaceae',
	
	#STEP2: CREATE DICT OF DICT {Kingdom:{'k_Bacteria':['Sample_1', 'Sample_2']}}
	TAX_metadata_DataFrame_DICT = {}
	TAX_metadata_DATAFrame_header_LIST = []
	TAX_metadata_DataFrame_DICT['TAX_ALIAS'] = []
	TAX_metadata_DataFrame_DICT['TAX_ID'] = []
	TAX_metadata_DATAFrame_header_LIST = ['TAX_ALIAS', 'TAX_ID']

	TAX_map_DICT = {}
	alias_counter = 1
	TAX_temp_DICT = {}
	for each_lineage in lineage_Name_LIST:
		TAX_map_DICT[each_lineage] = {}
		for each_taxonomy in lineage_target_DICT[each_lineage]:
			TAX_alias_name = 'TAX_ALIAS_' + str(alias_counter).zfill(6)
			TAX_metadata_DataFrame_DICT['TAX_ALIAS'].append(TAX_alias_name)
			TAX_metadata_DataFrame_DICT['TAX_ID'].append(each_taxonomy)
			
			TAX_temp_DICT[each_lineage + ';;' + each_taxonomy] = TAX_alias_name

			TAX_map_DICT[each_lineage][TAX_alias_name] = []
			alias_counter += 1
	#print TAX_metadata_DataFrame_DICT

	for each_OTU_ALIAS, each_taxonomy in itertools.izip(OTU_metadata_DataFrame_DICT['OTU_ALIAS'], OTU_metadata_DataFrame_DICT['taxonomy']):

		taxonomy_LIST = each_taxonomy.split(';')
		#print taxonomy_LIST
		########################################
		TAX_map_DICT['Kingdom'][TAX_temp_DICT['Kingdom;;' + taxonomy_LIST[0]]].append(each_OTU_ALIAS)
		#############################################################################################################
		TAX_map_DICT['Phylum'][TAX_temp_DICT['Phylum;;' + taxonomy_LIST[1]]].append(each_OTU_ALIAS)
		#############################################################################################################
		TAX_map_DICT['Class'][TAX_temp_DICT['Class;;' + taxonomy_LIST[2]]].append(each_OTU_ALIAS)
		#############################################################################################################
		TAX_map_DICT['Order'][TAX_temp_DICT['Order;;' + taxonomy_LIST[3]]].append(each_OTU_ALIAS)
		#############################################################################################################
		TAX_map_DICT['Family'][TAX_temp_DICT['Family;;' + taxonomy_LIST[4]]].append(each_OTU_ALIAS)
		#############################################################################################################
		TAX_map_DICT['Genus'][TAX_temp_DICT['Genus;;' + taxonomy_LIST[5]]].append(each_OTU_ALIAS)
		#############################################################################################################
		TAX_map_DICT['Species'][TAX_temp_DICT['Species;;' + taxonomy_LIST[6]]].append(each_OTU_ALIAS)
		#############################################################################################################
	#print TAX_map_DICT
	#{'Kingdom': {'TAX_ALIAS_000001': ['OTU_ALIAS_000001', 'OTU_ALIAS_000002',
	#'Class': {'TAX_ALIAS_000012': ['OTU_ALIAS_000032'], 'TAX_ALIAS_000013': ['OTU_ALIAS_000074'

	#STEP 3: MULTI_TAX_SHARED_DICT = {}
	OTU_name_LIST = []
	MULTI_TAX_SHARED_DICT = {}
	for each_lineage in lineage_Name_LIST:
		OTU_name_LIST = ['label', 'Groups', 'numOtus']
		MULTI_TAX_SHARED_DICT[each_lineage] = {}
		MULTI_TAX_SHARED_DICT[each_lineage]['label'] = SHARED_file_DataFrame_DICT['label']
		MULTI_TAX_SHARED_DICT[each_lineage]['Groups'] = SHARED_file_DataFrame_DICT['Groups']
		MULTI_TAX_SHARED_DICT[each_lineage]['numOtus'] = [len(TAX_map_DICT[each_lineage])] * len(SHARED_file_DataFrame_DICT['Groups'])
		OTU_name_LIST.extend(TAX_map_DICT[each_lineage])
		for each_TAX_alias in TAX_map_DICT[each_lineage]:
			MULTI_TAX_SHARED_DICT[each_lineage][each_TAX_alias] = []
			list_of_values = []
			for each_OTU_ALIAS in TAX_map_DICT[each_lineage][each_TAX_alias]:
					list_of_values.append(map(int, SHARED_file_DataFrame_DICT[each_OTU_ALIAS]))
			MULTI_TAX_SHARED_DICT[each_lineage][each_TAX_alias] = [sum(x) for x in zip(*list_of_values)]
		lineage_data_frame = pandas.DataFrame.from_dict(MULTI_TAX_SHARED_DICT[each_lineage])
		lineage_data_frame.to_excel(EXCEL_writer, sheet_name=each_lineage, columns=OTU_name_LIST, index=False, engine='xlsxwriter')
	#print MULTI_TAX_SHARED_DICT['Phylum']
	#{'TAX_ALIAS_000005': ['1', '1', '0', '0', '0', '5', '0', '3', '2'], 'TAX_ALIAS_000004': ['7', '2', '0', '0', '0', '0', '0', '3', '5'], 'TAX_ALIAS_000007': ['1'

	#Writing Metadata
	TAX_metadata_data_frame = pandas.DataFrame.from_dict(TAX_metadata_DataFrame_DICT)
	TAX_metadata_data_frame.to_excel(EXCEL_writer, sheet_name='TAX_metadata', columns=TAX_metadata_DATAFrame_header_LIST, index=False, engine='xlsxwriter')

	EXCEL_writer.save()
	return True


def match_two_list(list_A, list_B):

	return set(list_A) & set(list_B)


def excel_to_vertical_dict_converter(Excel_file_PATH, SHEETname):
	DATAFRAME = pandas.read_excel(Excel_file_PATH, sheetname=SHEETname)
	any_vertical_LIST = []
	any_vertical_LIST = DATAFRAME.columns.values.tolist()
	any_vertical_DICT = {}
	for each_column in any_vertical_LIST:
		each_column_value_list = DATAFRAME[each_column].values.tolist()
		any_vertical_DICT[str(each_column)] = map(str, each_column_value_list)
	any_vertical_LIST = map(str, any_vertical_LIST)
	return (any_vertical_DICT, any_vertical_LIST)


def excel_to_matrix_converter(Excel_file_PATH, SHEETname, index_col=None):
	if index_col is None:
		index_col = 0
	vertical_DICT, vertical_LIST = excel_to_vertical_dict_converter(Excel_file_PATH, SHEETname)
	Excel_MATRIX = {}
	Excel_LIST = vertical_LIST
	for each_SAMPLE in vertical_DICT[vertical_LIST[index_col]]:
		Excel_MATRIX[each_SAMPLE] = {}
		each_SAMPLE_index = vertical_DICT[vertical_LIST[index_col]].index(each_SAMPLE)
		for each_column in vertical_LIST:
			Excel_MATRIX[each_SAMPLE][each_column] = vertical_DICT[each_column][each_SAMPLE_index]
	return (Excel_MATRIX, Excel_LIST)


def interpret_excel_file(Excel_file_PATH, interpreted_Excel_file_PATH):
	EXCEL_writer = pandas.ExcelWriter(interpreted_Excel_file_PATH)
	TAXIC_Name_LIST = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

	TAX_metadata_MATRIX, TAX_metadata_LIST = excel_to_matrix_converter(Excel_file_PATH, 'TAX_metadata', 0)
	Sample_metadata_MATRIX, Sample_metadata_LIST = excel_to_matrix_converter(Excel_file_PATH, 'Sample_metadata', 0)
	for each_TAX_level in TAXIC_Name_LIST:
		TAXIC_data_vertical_DICT, TAXIC_data_vertical_LIST = excel_to_vertical_dict_converter(Excel_file_PATH, each_TAX_level)
		new_TAXIC_DICT = {}
		new_TAXIC_COLUMNS_LIST = []
		for each_TAX_column in TAXIC_data_vertical_LIST:
			#print each_TAX_column

			if each_TAX_column in TAX_metadata_MATRIX:
				#print TAX_metadata_MATRIX[each_TAX_column]
				new_TAXIC_DICT[TAX_metadata_MATRIX[each_TAX_column]['TAX_ID']] = TAXIC_data_vertical_DICT[each_TAX_column]
				new_TAXIC_COLUMNS_LIST.append(TAX_metadata_MATRIX[each_TAX_column]['TAX_ID'])
			else:
				new_TAXIC_DICT[each_TAX_column] = TAXIC_data_vertical_DICT[each_TAX_column]

				new_TAXIC_COLUMNS_LIST.append(each_TAX_column)

		NEW_TAXIC_data_frame = pandas.DataFrame.from_dict(new_TAXIC_DICT)
		NEW_TAXIC_data_frame.to_excel(EXCEL_writer, sheet_name=each_TAX_level, columns=new_TAXIC_COLUMNS_LIST, index=False, engine='xlsxwriter')
	EXCEL_writer.save()
	return True


def remove_extension_files(outputdir, extension):
	extension_list = []
	extension_list.append(extension)
	scanned_container = []
	flag = scandirs(outputdir, scanned_container, extension_list, 'partial')
	print "Scanning to find", extension, "Files.."
	if flag is False:
		print "Failed :("
		print "This extension is not available: ", extension_list
		sys.exit(2)
	else:
		print "NICE :) Found them"
		counter = 1
		for file in scanned_container:
			#print "Removing File#", str(counter), ":", file
			counter += 1
			check_it_and_remove_it(file)
	return True
# ################################# MAFFT FUNCTIONS ########################## #


def mafft_align(processors, outputdir, stderr, stdout, run_pid, mafft_exec_path, fasta_file, mafft_out, mafft_parameter):
	space = ' '
	mafft_string = 'nohup' + space
	mafft_string += mafft_exec_path + space
	mafft_string += mafft_parameter + '--thread' + space + processors + space
	mafft_string += fasta_file + space
	#mafft_string += '> ' + mafft_out + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	mafft_string += '> ' + mafft_out + ' & echo $! > ' + run_pid
	exec_dict = {}
	print mafft_string
	report(mafft_string)
	exec_dict = execute([mafft_string])
	if exec_dict["exitCode"] != 0:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	else:
		return True


def test_mafft(processors, outputdir, stderr, stdout, run_pid, mafft_exec_path):
	space = ' '
	mafft_string = 'nohup' + space + mafft_exec_path + space + '--help' + space + '--thread' + space + processors + space + '> ' + stderr + ' 2> ' + stdout + ' & echo $! > ' + run_pid
	print "EXECUTING: ", mafft_string
	report(mafft_string)
	exec_dict = {}
	exec_dict = execute([mafft_string])
	if exec_dict["exitCode"] != 0:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	else:
		return True


def process_taxonomy_file(taxonomy_file_PATH, processed_taxonomy_file_PATH):
	f = open(taxonomy_file_PATH, 'rU')
	new_tax_string = ''
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		TAX_ID = line[0]
		TAX_taxonomy = line[1]
		TAX_taxonomy = remove_ambiguity_from_taxonomy(TAX_taxonomy.split(';'))
		TAX_taxonomy = TAX_taxonomy[:-1] + '_TAG_' + TAX_ID + ';'
		new_tax_string += TAX_ID + '\t' + TAX_taxonomy + '\n'
	f.close()
	write_string_down(new_tax_string, processed_taxonomy_file_PATH)
	return True


def remove_ambiguity_from_taxonomy(taxonomy_LIST):
	taxonomy_identifier = ['k_', 'p_', 'c_', 'o_', 'f_', 'g_', 's_']
	#Step1:Clean up
	taxonomy_LIST = filter(None, taxonomy_LIST)
	Clean_taxonomy_List = []
	for each_tax_level in taxonomy_LIST:
		each_tax_level = str(each_tax_level).lower()
		each_tax_level = slugify(each_tax_level)
		if each_tax_level in ['root', 'unknown', 'domain', 'k_', 'p_', 'c_', 'o_', 'f_', 'g_', 's_']:
			continue
		if each_tax_level[:2] in taxonomy_identifier:
			each_tax_level = each_tax_level[2:]
		each_tax_level = each_tax_level.title()
		Clean_taxonomy_List.append(each_tax_level)
	#Step2: Length_Trimming
	if len(Clean_taxonomy_List) > 7:
		Clean_taxonomy_List = Clean_taxonomy_List[:7]

	last_index = 0
	tax_string = ''
	Last_tax = ''

	for each_level in taxonomy_identifier:
		each_level_index = taxonomy_identifier.index(each_level)
		try:
			Clean_taxonomy_List[each_level_index]
			last_index = each_level_index
			last_level = each_level
			if Clean_taxonomy_List[last_index][0:2] != each_level:
				tax_string += last_level + Clean_taxonomy_List[last_index] + ';'
				Last_tax = last_level + Clean_taxonomy_List[last_index] + ';'
			else:
				tax_string += Clean_taxonomy_List[last_index] + ';'
				Last_tax = Clean_taxonomy_List[last_index] + ';'
		except IndexError:
			tax_string += '*' + Last_tax
	return tax_string.replace('_;', ';')
# ################################### MOTHUR FUNCTIONS ####################### #


def test_mothur(processors, outputdir, stderr, stdout, run_pid, mothur_exec_PATH):
	# test mothur to see if it is working and grab the version of mothur by scanning its output log
	mothur_input_dictionary = {}
	command = 'get.current'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_PATH
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	#parameter_list = []
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_make_file(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, input_file_directory):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	space = ' '
	command = 'make.file'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('inputdir=' + input_file_directory)
	parameter_list.append(',' + space + 'prefix=miseq_sop')
	mothur_input_dictionary['parameters'] = parameter_list
	#here we construct the object
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_make_contigs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, files_PATH):
	mothur_input_dictionary = {}
	#space = ' '
	command = 'make.contigs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('file=' + files_PATH)
	mothur_input_dictionary['parameters'] = parameter_list
	#here we construct the object
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_screen_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, group_file, maxambig, maxhomop, min_len):
	mothur_input_dictionary = {}
	command = 'screen.seqs'
	space = ' '
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'group=' + group_file)
	if maxambig is not False:
		parameter_list.append(',' + space + 'maxn=' + maxambig)
	if min_len is not False:
		parameter_list.append(',' + space + 'minlength=' + min_len)
		parameter_list.append(',' + space + 'maxlength=' + min_len)
	if maxhomop is not False:
		parameter_list.append(',' + space + 'maxhomop=' + maxhomop)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_summary_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	command = 'summary.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_unique_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	command = 'unique.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_count_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, group_file, name_file):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	space = ' '
	command = 'count.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('group=' + group_file)
	parameter_list.append(',' + space + 'name=' + name_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_precluster(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, count_file):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	space = ' '
	command = 'pre.cluster'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'count=' + count_file)
	parameter_list.append(',' + space + 'diffs=2')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_chimera_vsearch(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, count_file):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	space = ' '
	command = 'chimera.vsearch'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'count=' + count_file)
	parameter_list.append(',' + space + 'dereplicate=t')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_remove_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, count_file, accnos_file):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	space = ' '
	command = 'remove.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'count=' + count_file)
	parameter_list.append(',' + space + 'accnos=' + accnos_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_classify_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, reference, taxonomy):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'classify.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'reference=' + reference)
	parameter_list.append(',' + space + 'taxonomy=' + taxonomy)
	parameter_list.append(',' + space + 'probs=F, method=knn, cutoff=80, iters=100, numwanted=1')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_dist_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, cutoff_limit):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'dist.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'cutoff=' + cutoff_limit)
	parameter_list.append(',' + space + 'output=column')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_cluster_split(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, distance_file, count_table):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'cluster.split'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('column=' + distance_file)
	parameter_list.append(',' + space + 'count=' + count_table + ',' + space + 'large=T, method=average, cutoff=0.03')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_make_shared(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, list_file, count_table):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'make.shared'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('list=' + list_file)
	parameter_list.append(',' + space + 'count=' + count_table + ',' + space + 'label=0.03')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_remove_rare(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, list_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'remove.rare'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('list=' + list_file)
	parameter_list.append(',' + space + 'nseqs=2')
	parameter_list.append(',' + space + 'label=0.03')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_get_oturep(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, list_file, count_table, distance_file, fasta_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'get.oturep'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('column=' + distance_file)
	parameter_list.append(',' + space + 'list=' + list_file)
	parameter_list.append(',' + space + 'fasta=' + fasta_file)
	parameter_list.append(',' + space + 'count=' + count_table)
	parameter_list.append(',' + space + 'label=0.03, large=true, method=distance')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_align_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, reference_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'align.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('candidate=' + fasta_file)
	parameter_list.append(',' + space + 'template=' + reference_file)
	parameter_list.append(',' + space + 'flip=t')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_get_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, accnos_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'get.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'accnos=' + accnos_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True
# ################################### MAIN FUNCTION ########################## #


def main(argv):

	report_string = ''
	# ++++++++++++++++++++++++++++++ PARSE INPUT ARGUMENTS
	parser = argparse.ArgumentParser()
	main_file = parser.add_argument_group('Main file parameters')
	main_file.add_argument("--input", help="zipped folder contains fastq files", action='store')
	args = parser.parse_args()

	# ++++++++++++++++++++++++++++++ BEURACRATICS PROCEDURES
	report_string += "######################################################################################################################################\n"
	print "######################################################################################################################################"
	report_string += "MISEQ SOP 1.0 EXECUTION HAS INITIATED" + '\n'
	print "MISEQ SOP 1.0 EXECUTION HAS INITIATED"
	report_string += "Initiation time: " + time.strftime("%Y-%m-%d %H:%M:%S") + '\n'
	print "Initiation time: ", time.strftime("%Y-%m-%d %H:%M:%S")
	report_string += "###################################################################" + '\n'
	print "###################################################################"
	report_string += "INFORMATION ABOUT THE ENVIRONMENT, EXECUTABLES AND PROVIDED DATA" + '\n'
	print "INFORMATION ABOUT THE ENVIRONMENT, EXECUTABLES AND PROVIDED DATA"
	report_string += "###################################################################" + '\n'
	print "###################################################################"
	report_string += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" + '\n'
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "COMMAND LINE:"
	report_string += "COMMAND LINE:\n"
	report_string += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" + '\n'
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	commandline_string = 'python ' + ' '.join(sys.argv) + '\n'
	print commandline_string
	report_string += commandline_string + '\n'
	report_string += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" + '\n'
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report_string += "ARGUMENTS:" + '\n'
	print "ARGUMENTS:"
	report_string += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" + '\n'
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

	# ++++++++++++++++++++++++++++++ OUTPUT DIRECTORY CHECKING
	args.outputdir = DEFAULT_OUTPUTDIR
	global report_file
	report_file = args.outputdir + "miseq_sop_report.txt"
	check_it_and_remove_it(report_file, True)
	report(report_string)
	# ------------------------------ END OF OUTPUT DIRECTORY CHECKING

	# ++++++++++++++++++++++++++++++ PROCESSORS CHECKING
	args.processors = DEFAULT_PROCESSORS
	# ------------------------------ END OF PROCESSORS CHECKING

	# ++++++++++++++++++++++++++++++ PREFIX NAME CHECKING
	args.prefix = DEFAULT_PREFIX
	# ------------------------------ END OF PREFIX NAME CHECKING

	# ++++++++++++++++++++++++++++++ EXECUTIVE DIRECTORY CHECKING
	args.execdir = DEFAULT_EXECDIR
	# ------------------------------ END OF EXECUTIVE DIRECTORY CHECKING

	# ++++++++++++++++++++++++++++++ TAXONOMY DIRECTORY CHECKING
	args.taxdir = DEFAULT_TAXDIR
	# ------------------------------ END OF TAXONOMY DIRECTORY CHECKING

	# ++++++++++++++++++++++++++++++ CHECKING EXECUTIVES
	print "\n###################################################################"
	report("\n###################################################################")
	print "VERIFYING THE SANITY/VERSION OF EXECUTABLES IN EXECUTIVE DIRECTORY[--execdir]"
	report("VERIFYING THE SANITY/VERSION OF EXECUTABLES IN EXECUTIVE DIRECTORY[--execdir]")
	print "###################################################################\n"
	report("###################################################################\n")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("0: ENVIRONMENT")
	print "0: ENVIRONMENT"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "Operating System version is: ", platform.platform()
	report("Operating System version is: " + platform.platform())
	python_version = sys.version.split(' (')[0]
	print "Python version is: ", python_version
	report("Python version is: " + python_version)
	if float(python_version[0:3]) < 2.7:
		error("Python version is older than 2.7")
		print "Python version is older than 2.7"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("1: MOTHUR")
	print "1: MOTHUR"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	mothur_exec_PATH = args.execdir + 'mothur/mothur'
	if isFileExist(mothur_exec_PATH) is False:
		error("Your mothur file path has Access/Exist issue")
		print "Your mothur file path has Access/Exist issue"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		print "mothur execution file path is: ", mothur_exec_PATH
		report("mothur execution file path is: " + mothur_exec_PATH)
		#report("Testing mothur executables: ")
		print "Testing mothur executable: "
		report("Testing mothur executable: ")
		flag, stderr = execute_functions(test_mothur, args.processors, args.outputdir, 'multi', 'mothur', mothur_exec_PATH)
		if flag is False:
			error("[" + FAILED_MARK + "]")
			print "[", FAILED_MARK, "]"
			print "Execution of mothur failed!!!"
			error("Execution of mothur failed!!!")
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			target_lines = []
			flag = parse_mothur_logfile(stderr, 'mothur v.', target_lines)
			if flag is False:
				print "This keyword is not avalaible: mothur v."
				error("This keyword is not avalaible: mothur v.")
				print "ABORTING!!!"
				error("ABORTING!!!")
				sys.exit(2)
			report("Mothur executables responding successfully!!!")
			print "Mothur executables responding successfully!!!"
			report("version of mothur executables: " + target_lines[0])
			print "version of mothur executables:", target_lines[0]
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("2: MAFFT")
	print "2: MAFFT"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	mafft_exec_path = args.execdir + 'mafft.bat'
	if isFileExist(mafft_exec_path) is False:
		error("Your mafft path has Access/Exist issue")
		print "Your mafft path has Access/Exist issue"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("mafft execution file is: " + mafft_exec_path)
		print "mafft execution file is: ", mafft_exec_path
		report("Testing mafft executables: ")
		print "Testing mafft executables:"
		flag, stderr = execute_functions(test_mafft, args.processors, args.outputdir, 'multi', 'mafft', mafft_exec_path)
		if flag is False:
			print "Execution of mafft failed!!!"
			error("Execution of mafft failed!!!")
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			report("mafft executables responding successfully!!!")
			print "mafft executables responding successfully!!!"
			report("version of mafft executables: " + stderr.split('\n')[1])
			print "version of mafft executables:", stderr.split('\n')[1]
	# ------------------------------ END OF CHECKING EXECUTIVES
	# ------------------------------ END OF BEURACRATICS PROCEDURES
	
	# ++++++++++++++++++++++++++++++ CHECKING INPUTS
	print "\n###################################################################"
	report("\n###################################################################")
	print "VERIFYING THE SANITY/VALIDITY OF INPUT FILES"
	report("VERIFYING THE SANITY/VALIDITY OF INPUT FILES")
	print "###################################################################\n"
	report("###################################################################\n")
	print "EXTRACING INPUT is in progress"
	report("EXTRACING INPUT is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	if args.input is not None:
		flag = zip_file_validation_process(args.input, DEFAULT_INPUTDIR)
		if flag is True:
			print "EXTRACING INPUT STEP: PASSED!!!"
			report("EXTRACING INPUT STEP: PASSED!!!")
			args.input = DEFAULT_INPUTDIR
		else:
			error("EXTRACING INPUT STEP: FAILED!!!")
			print "EXTRACING INPUT STEP: FAILED!!!"
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
	else:
		args.input = DEFAULT_TESTDIR
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	
	# #####################################################################
	# END OF BEURACRATICS PROCEDURES
	# #####################################################################
	
	# #####################################################################
	# INPUT DATA PREPARATION
	# #####################################################################
	# +++++++++++++++++++++++++++++ INPUT DATA PREPARATION
	print "INPUT DATA PREPARATION is in progress"
	report("INPUT DATA PREPARATION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	files_path = args.outputdir + 'miseq_sop_files_path_STEP1.txt'
	flag = data_prep(args.input, files_path, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "INPUT DATA PREPARATION STEP: PASSED!!!"
		report("INPUT DATA PREPARATION STEP: PASSED!!!")
	else:
		error("INPUT DATA PREPARATION STEP: FAILED!!!")
		print "INPUT DATA PREPARATION STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF INPUT DATA PREPARATION

	# +++++++++++++++++++++++++++++ MAKE CONTIGS
	print "MAKE CONTIGS is in progress"
	report("MAKE CONTIGS is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	contigs_fasta = args.outputdir + 'miseq_sop_contigs_fasta_STEP2.fasta'
	contigs_groups = args.outputdir + 'miseq_sop_contigs_groups_STEP2.txt'
	flag = merge_relabel(files_path, contigs_fasta, contigs_groups, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "MAKE CONTIGS STEP: PASSED!!!"
		report("MAKE CONTIGS STEP: PASSED!!!")
	else:
		error("MAKE CONTIGS STEP: FAILED!!!")
		print "MAKE CONTIGS STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF MAKE CONTIGS

	# ###########################################################
	# DATA FILTERING
	# ###########################################################

	# +++++++++++++++++++++++++++++ FILTER FASTA FILES
	print "FILTER FASTA FILES is in progress"
	report("FILTER FASTA FILES is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	filtered_fasta_file = args.outputdir + 'miseq_sop_filtered_fasta_STEP3.fasta'
	filtered_group_file = args.outputdir + 'miseq_sop_filtered_group_STEP3.fasta'
	flag = filter_fasta_file(contigs_fasta, contigs_groups, filtered_fasta_file, filtered_group_file, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "FILTER FASTA FILES STEP: PASSED!!!"
		report("FILTER FASTA FILES STEP: PASSED!!!")
	else:
		error("FILTER FASTA FILES STEP: FAILED!!!")
		print "FILTER FASTA FILES STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF FILTER FASTA FILES

	# +++++++++++++++++++++++++++++ REMOVING DUPLICATES
	print "REMOVING DUPLICATES is in progress"
	report("REMOVING DUPLICATES is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	unique_fasta_file = args.outputdir + 'miseq_sop_unique_fasta_STEP4.fasta'
	unique_names_file = args.outputdir + 'miseq_sop_unique_names_STEP4.txt'
	flag = removing_duplicates(filtered_fasta_file, unique_fasta_file, unique_names_file, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "REMOVING DUPLICATES STEP: PASSED!!!"
		report("REMOVING DUPLICATES STEP: PASSED!!!")
	else:
		error("REMOVING DUPLICATES STEP: FAILED!!!")
		print "REMOVING DUPLICATES STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF REMOVING DUPLICATES

	# +++++++++++++++++++++++++++++ COUNT SEQS FILE
	print "COUNT SEQS FILE is in progress"
	report("COUNT SEQS FILE is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	count_seqs_file = args.outputdir + 'miseq_sop_count_seqs_STEP5.txt'
	flag = counting_sequences(filtered_group_file, unique_names_file, count_seqs_file, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "COUNT SEQS FILE STEP: PASSED!!!"
		report("COUNT SEQS FILE STEP: PASSED!!!")
	else:
		error("COUNT SEQS FILE STEP: FAILED!!!")
		print "COUNT SEQS FILE STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF COUNT SEQS FILE

	# +++++++++++++++++++++++++++++ PRECLUSTER
	print "PRECLUSTER is in progress"
	report("PRECLUSTER is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	precluster_fasta = args.outputdir + 'miseq_sop_precluster_fasta_STEP6.fasta'
	precluster_count = args.outputdir + 'miseq_sop_precluster_count_STEP6.txt'
	flag = precluster(unique_fasta_file, count_seqs_file, precluster_fasta, precluster_count, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "PRECLUSTER STEP: PASSED!!!"
		report("PRECLUSTER STEP: PASSED!!!")
	else:
		error("PRECLUSTER STEP: FAILED!!!")
		print "PRECLUSTER STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF PRECLUSTER

	# +++++++++++++++++++++++++++++ CHIMERA VSEARCH
	print "CHIMERA VSEARCH is in progress"
	report("CHIMERA VSEARCH is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	chimera_report_file = args.outputdir + 'miseq_sop_chimeric_report_file_STEP7.txt'
	flag = chimera_scanning(precluster_fasta, precluster_count, chimera_report_file, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "CHIMERA VSEARCH STEP: PASSED!!!"
		report("CHIMERA VSEARCH STEP: PASSED!!!")
	else:
		error("CHIMERA VSEARCH STEP: FAILED!!!")
		print "CHIMERA VSEARCH STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF CHIMERA VSEARCH

	# +++++++++++++++++++++++++++++ REMOVING CHIMEREIC SEQUENCES
	print "REMOVING CHIMEREIC SEQUENCES is in progress"
	report("REMOVING CHIMEREIC SEQUENCES is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	chimeric_removed_fasta = args.outputdir + 'miseq_sop_chimeric_removed_fasta_STEP8.fasta'
	chimeric_removed_count = args.outputdir + 'miseq_sop_chimeric_removed_count_STEP8.txt'
	flag = remove_chimeric_seqs(precluster_fasta, precluster_count, chimera_report_file, chimeric_removed_fasta, chimeric_removed_count, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "REMOVING CHIMEREIC SEQUENCES STEP: PASSED!!!"
		report("REMOVING CHIMEREIC SEQUENCES STEP: PASSED!!!")
	else:
		error("REMOVING CHIMEREIC SEQUENCES STEP: FAILED!!!")
		print "REMOVING CHIMEREIC SEQUENCES STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF REMOVING CHIMEREIC SEQUENCES

	# ########################################################
	# UPGMA CLUSTERING
	# ########################################################

	# +++++++++++++++++++++++++++++ MULTIPLE ALIGNMENT USING MAFFT
	print "MULTIPLE ALIGNMENT USING MAFFT is in progress"
	report("MULTIPLE ALIGNMENT USING MAFFT is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	space = ' '
	mafft_parameter = ''
	mafft_parameter = '--parttree' + space
	multiple_aligned_fasta = args.outputdir + 'miseq_sop_mafft_multiple_aligned_STEP9.fasta'

	flag, stderr = execute_functions(mafft_align, args.processors, args.outputdir, 'multi', 'mafft', mafft_exec_path, chimeric_removed_fasta, multiple_aligned_fasta, mafft_parameter)
	if flag is True:
		print "MULTIPLE ALIGNMENT USING MAFFT STEP: PASSED!!!"
		report("MULTIPLE ALIGNMENT USING MAFFT STEP: PASSED!!!")
	else:
		error("MULTIPLE ALIGNMENT USING MAFFT STEP: FAILED!!!")
		print "MULTIPLE ALIGNMENT USING MAFFT STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF MULTIPLE ALIGNMENT USING MAFFT

	# +++++++++++++++++++++++++++++ DISTANCE MATRIX CALCULATION
	print "DISTANCE MATRIX CALCULATION is in progress"
	report("DISTANCE MATRIX CALCULATION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	distance_matrix_file = args.outputdir + 'miseq_sop_distance_matrix_STEP10.txt'
	flag = distance_matrix(multiple_aligned_fasta, distance_matrix_file, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "DISTANCE MATRIX CALCULATION STEP: PASSED!!!"
		report("DISTANCE MATRIX CALCULATION STEP: PASSED!!!")
	else:
		error("DISTANCE MATRIX CALCULATION STEP: FAILED!!!")
		print "DISTANCE MATRIX CALCULATION STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	report("####################################################################")
	print "####################################################################"
	# ----------------------------- END OF DISTANCE MATRIX CALCULATION

	# +++++++++++++++++++++++++++++ UPGMA CLUSTERING
	print "UPGMA CLUSTERING is in progress"
	report("UPGMA CLUSTERING is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	cluster_list_file = args.outputdir + 'miseq_sop_cluster_list_file_STEP11.txt'
	flag = cluster_split(distance_matrix_file, chimeric_removed_count, cluster_list_file, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "UPGMA CLUSTERING STEP: PASSED!!!"
		report("UPGMA CLUSTERING STEP: PASSED!!!")
	else:
		error("UPGMA CLUSTERING STEP: FAILED!!!")
		print "UPGMA CLUSTERING STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	report("####################################################################")
	print "####################################################################"
	# ----------------------------- END OF UPGMA CLUSTERING

	# +++++++++++++++++++++++++++++ REMOVING LOW ABUNDANCE OTUS
	print "REMOVING LOW ABUNDANCE OTUS is in progress"
	report("REMOVING LOW ABUNDANCE OTUS is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	filtered_list_file = args.outputdir + 'miseq_sop_filtered_list_file_STEP12.txt'
	flag = filter_list_file(cluster_list_file, filtered_list_file, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "MAKING SHARED FILE STEP: PASSED!!!"
		report("MAKING SHARED FILE STEP: PASSED!!!")
	else:
		error("MAKING SHARED FILE STEP: FAILED!!!")
		print "MAKING SHARED FILE STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	report("####################################################################")
	print "####################################################################"
	# ----------------------------- END OF REMOVING LOW ABUNDANCE OTUS

	# +++++++++++++++++++++++++++++ OTU REPRESNTATIVE SEQUENCE
	print "OTU REPRESNTATIVE SEQUENCE is in progress"
	report("OTU REPRESNTATIVE SEQUENCE is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	representative_OTU_fasta_file_PATH = args.outputdir + 'miseq_sop_representative_OTU_fasta_file_STEP13.fasta'
	representative_OTU_count_file_PATH = args.outputdir + 'miseq_sop_representative_OTU_count_file_STEP13.txt'
	flag = get_representative_OTU(filtered_list_file, chimeric_removed_count, distance_matrix_file, chimeric_removed_fasta, representative_OTU_fasta_file_PATH, representative_OTU_count_file_PATH, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "OTU REPRESNTATIVE SEQUENCE STEP: PASSED!!!"
		report("OTU REPRESNTATIVE SEQUENCE STEP: PASSED!!!")
	else:
		error("OTU REPRESNTATIVE SEQUENCE STEP: FAILED!!!")
		print "OTU REPRESNTATIVE SEQUENCE STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	report("####################################################################")
	print "####################################################################"
	# ----------------------------- END OF OTU REPRESNTATIVE SEQUENCE

	# +++++++++++++++++++++++++++++ MAKING CENTROIDS FASTA
	print "MAKING CENTROIDS FASTA is in progress"
	report("MAKING CENTROIDS FASTA is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	centroids_fasta = args.outputdir + 'miseq_sop_UPGMA_centroids_STEP14.fasta'
	check_it_and_remove_it(centroids_fasta)
	tag = 'OTU'
	suffix = ''
	flag = relable_fasta(representative_OTU_fasta_file_PATH, tag, suffix, centroids_fasta)
	if flag is False:
		print "Execution of relabel fasta failed!!!"
		sys.exit(2)
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	report("####################################################################")
	print "####################################################################"
	# ----------------------------- END OF MAKING CENTROIDS FASTA

	# ########################################################
	# CLASSIFICATION OF CENTROIDS
	# ########################################################
	
	# +++++++++++++++++++++++++++++ CLASSIFICATION OF CENTROIDS
	print "CLASSIFICATION OF CENTROIDS is in progress"
	report("CLASSIFICATION OF CENTROIDS is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	centroids_classification_report = args.outputdir + 'miseq_sop_centroids_classification_report_STEP15.txt'
	centroids_classification_fasta = args.outputdir + 'miseq_sop_centroids_classification_fasta_STEP15.fasta'
	RDP_reference = args.taxdir + 'gg_13_8_99.fasta'
	RDP_taxonomy = args.taxdir + 'gg_13_8_99.gg.tax'
	flag = classify_centroids(centroids_fasta, RDP_reference, RDP_taxonomy, centroids_classification_report, centroids_classification_fasta, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "CLASSIFICATION OF CENTROIDS STEP: PASSED!!!"
		report("CLASSIFICATION OF CENTROIDS STEP: PASSED!!!")
	else:
		error("CLASSIFICATION OF CENTROIDS STEP: FAILED!!!")
		print "CLASSIFICATION OF CENTROIDS STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)

	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	report("####################################################################")
	print "####################################################################"
	# ----------------------------- END OF CLASSIFICATION OF CENTROIDS

	# +++++++++++++++++++++++++++++ GLOBAL ALIGNMENT
	print "GLOBAL ALIGNMENT is in progress"
	report("GLOBAL ALIGNMENT is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	global_alignment_report = args.outputdir + 'miseq_sop_centroids_global_alignment_report_STEP16.txt'
	flag = global_alignment(contigs_fasta, centroids_classification_fasta, global_alignment_report, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		print "GLOBAL ALIGNMENT STEP: PASSED!!!"
		report("GLOBAL ALIGNMENT STEP: PASSED!!!")
	else:
		error("GLOBAL ALIGNMENT STEP: FAILED!!!")
		print "GLOBAL ALIGNMENT STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)

	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	report("####################################################################")
	print "####################################################################"
	# ----------------------------- END OF GLOBAL ALIGNMENT
	
	# +++++++++++++++++++++++++++++ GENERATE RESULT
	print "GENERATE RESULT is in progress"
	report("GENERATE RESULT is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	#shared_file_PATH = args.outputdir + 'miseq_sop_result.tsv'
	shared_file_PATH = args.outputdir + 'miseq_sop_result_mothur_format.txt'
	SAS_file_PATH = args.outputdir + 'miseq_sop_result_sas_format.txt'
	BIOM_file_PATH = args.outputdir + 'miseq_sop_result_biom_format.biom'
	Excel_file_PATH = args.outputdir + 'miseq_sop_result_excel_format.xlsx'
	alignment_report_parser(global_alignment_report, centroids_classification_report, contigs_groups, Excel_file_PATH, SAS_file_PATH, BIOM_file_PATH, shared_file_PATH, args.outputdir)
	copy_file(BIOM_file_PATH, CURRENT_PATH + 'miseq_sop_result_biom_format.biom')
	copy_file(SAS_file_PATH, CURRENT_PATH + 'miseq_sop_result_sas_format.txt')
	copy_file(shared_file_PATH, CURRENT_PATH + 'miseq_sop_result_mothur_format.txt')

	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	report("####################################################################")
	print "####################################################################"
	# ----------------------------- END OF GENERATE RESULT
	flag = remove_extension_files(CURRENT_PATH, '.logfile')


	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "MISEQ SOP EXECUTION COMPLETED."
	report("MISEQ SOP EXECUTION COMPLETED.")
	report("Completion time: " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "Completion time: ", time.strftime("%Y-%m-%d %H:%M:%S")
# ################################### FINITO ################################# #
if __name__ == "__main__": main(sys.argv[1:])