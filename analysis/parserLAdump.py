#!/usr/bin/python
#
# Reformats the output of 'LAdump -cl' (provided in a file specified via argument 1)
# in bella format with one exception: read ID's are output in place of read names.
# The unique set of read IDs are recorded and output for lookup in the corresponding DAZZ_DB database instance.
# See translateDalignerOut.sh for the full translation procedure (which includes replacing the ID's with names).
#
# Uses a multi-reader single-writer parallelism strategy with a read-task pool and write-task queue. 
# The default number of threads is equal to the results of mp.cpu_count(). 
# Alternative values can be specified with the second argument.
#
# @author:mme
#
from sys import argv
import multiprocessing as mp

'''
 A pure function for parsing 3 'LAdump -cl' lines into the single corresponding BELLA line.

 Input: an iterable 'P', 'L', 'C' triplet
 Output: the input lines reformatted into a single BELLA alignment line
'''
def parse_las(lines):
  assert len(lines) == 3, "number of lines is %r instead of 3" % len(lines)
  assert lines[0].split()[0] == 'P', "bad chunking: first char of first line is not P"
  line=[]
  for l in lines:
    lsplit = l.split()
    if (lsplit[0]=='P'):
      line+=lsplit[1:4]
    if (lsplit[0]=='L'):
      line+=lsplit[1:]
    if (lsplit[0]=='C'):
      # r1 r2 c/n [3] [4] r1L [6] [7] r2L
      line.insert(3, lsplit[1])
      line.insert(4, lsplit[2])
      line.insert(6, lsplit[3])
      line.insert(7, lsplit[4])
  return line

'''
An asynchronous reader's task.

  infile; the file to read
  position; the position in the file from which to start reading
  size; the total number of liens to read
  queue; the write-queue to put formatted lines

  Returns: None
'''
def parse_las_chunk(infile, position, size, queue):
  with open(infile, 'r') as f:
    f.seek(position)
    sofar = 0
    lines = []
    while sofar < size:
      for _ in ['P','L','C']:
        lines.append(f.readline())
      parsed = parse_las(lines)
      queue.put( parsed ) # TODO submit lines in batches to reduce queue contention
      lines.clear()
      sofar += 3

'''
  The writer task.

  outfile; the name of the file to output reformatted lines to
  idfile; the name of the file to output read IDs 
  queue; the queue to get write tasks from
  finished; a message indicating all tasks are done (there will be nothing left to write)

  Return: None
'''
def write_lines(outfile, idfile, queue, finished):
  idset=set()
  with open(outfile, 'w') as out:
    while True:
      text = queue.get()
      if (text == finished): break
      idset.add(text[0])
      idset.add(text[1])
      formatted_line = ' '.join(x for x in text)
      out.write(formatted_line + '\n')
      out.flush() #TODO may be inefficient to flush after every line
# write the set of ids, 1 per line, to a file to use with DBshow
  with open(idfile, 'w') as out:
    for i in idset:
      out.write(i + '\n')
      out.flush() 

if __name__ == '__main__':
  if len(argv) < 2:
    print("1 argument is required: a file containing 'LAdump -cl' formatted alignment output")
    quit()
  
  threads = mp.cpu_count()
  if len(argv) == 3:
    threads = int(argv[2])

  manager = mp.Manager()
  pool = mp.Pool(threads+2)
  ln_queue = manager.Queue()

  # initialize reformatted-line writer
  finished = 'fin'
  ofilename=argv[1]+".bella"
  idfilename=argv[1]+".ids"
  writer = pool.apply_async(write_lines, (ofilename, idfilename, ln_queue, finished))

  #initialize readers/parsers
  tasks = []
  fname = argv[1]
  with open(fname, 'r') as f:
    num_algns = int(((f.readline()).split())[2]) # the total number of alignments is listed on the first line
    chunk_size = (num_algns // threads)*3
    last_chunk = chunk_size + ((num_algns % threads)*3)
    f.readline() # the first two lines of the file are header lines
    print("Num alignments, chunk_size, last_chunk")
    print(num_algns, chunk_size, last_chunk)
    if (chunk_size > 0):
      for i in range(threads-1):
        position = f.tell()
        tasks.append( pool.apply_async( parse_las_chunk(fname, position, chunk_size, ln_queue) ) )
        for _ in range(chunk_size): f.readline()
    if (last_chunk > 0):
      tasks.append( pool.apply_async( parse_las_chunk(fname, f.tell(), last_chunk, ln_queue)))

  # wait for all tasks to finish
  for task in tasks: 
    try:
      task.get()
    except TypeError as e:
      print("encountered known error invoking get() on finished task, continuing")

  ln_queue.put(finished)
  writer.get()

  # clean-up
  pool.close()

