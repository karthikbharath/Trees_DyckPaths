#!/usr/bin/env python

import os
import sys
import re
from optparse import OptionParser

#labels=[0,0,1,1,1,1,0,1,0,1,1,1,0,1,0,1,0,0,1,1,0,0,1,0,1,1,1,0,1,0,1,0,0,
 #       0,1,1,1,1,1,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1,1,0,0,0,1,1,1,1,1,0,0,0,
  #      1,1,1,1,1,0,0,1]

## These are fake label obtained from the one above through random permutation
labels=[0,0,0,1,1,1,0,1,1,0,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,1,1,0,1,0,1,0,0,
        0,1,01,0,0,0,1,1,1,0,0,0,1,1,1,1,0,1,0,1,1,0,0,0,1,0,1,1,1,0,0,0,
        0,1,0,1,1,0,0,1]


def parseArguments():
  parser = OptionParser('Process arguments')
  parser.add_option('-i', '--in-dir', 
                    dest='in_dir', 
                    default='.', 
                    action='store',
                    type='string',
                    help='input directory')
  parser.add_option('-o', '--out-dir', 
                    dest='out_dir', 
                    default='.', 
                    action='store',
                    type='string',
                    help='output directory')
  parser.add_option('-p', '--prefix', 
                    dest='prefix', 
                    default='result', 
                    action='store',
                    type='string',
                    help='prefix for relevant files')
  parser.add_option('-s', '--suffix', 
                    dest='suffix', 
                    default='txt', 
                    action='store',
                    type='string',
                    help='prefix for relevant files')
  (options, args) = parser.parse_args()

  return options

def main():
  options = parseArguments()

  if not os.path.isdir(options.in_dir): 
    raise InputError('Input directory invalid')

  if not os.path.isdir(options.out_dir):
    print '%s does not exist, creating it ... ' % options.out_dir,
    os.makedirs(options.out_dir)
    print 'done'

  for fname in os.listdir(options.in_dir):
    if fname.startswith(options.prefix) and fname.endswith(options.suffix):
      print 'Splitting %s ... ' % fname,
      fin     = open(options.in_dir + '/' + fname, 'r')
      f0T1out = open(options.out_dir + '/0-T1-' + fname, 'w')
      f1T1out = open(options.out_dir + '/1-T1-' + fname, 'w')
      f0T2out = open(options.out_dir + '/0-T2-' + fname, 'w')
      f1T2out = open(options.out_dir + '/1-T2-' + fname, 'w')
 
      for line in fin.readlines():
        # Skip the first line which only has header information
        if line.startswith('dist-type'): continue

        wholeTree = re.split(r' ', line)[2]
        treeCategory = re.split(r'_', wholeTree)[0]
        treeNumber = re.split(r'\.', re.split(r'_', wholeTree)[1])[0]
        treeLabel = labels[int(treeNumber)-1]

        if treeCategory == 'T1': 
          if 0 == treeLabel: print >>f0T1out, '%s' % line,
          else: print >>f1T1out, '%s' % line,
        else:
          if 0 == treeLabel: print >>f0T2out, '%s' % line,
          else: print >>f1T2out, '%s' % line,
 
      fin.close()
      f0T1out.close()
      f1T1out.close()
      f0T2out.close()
      f1T2out.close()
      print 'done'

# Main starts here
if __name__ == '__main__': main()
