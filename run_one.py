#!/usr/bin/env python
#
# Author: Anju Kambadur
# Please run as follows:
# ./run_one.py -n 1,2,3 4 -d <distribution> ...
#
# For help please see:
# ./run_one.py -h 
#
import os
import sys
import re
from optparse import OptionParser

def nodes_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, re.split(',', value))

def parseArguments():
  parser = OptionParser('Process arguments')
  parser.add_option('-d', '--distribution',
                    dest='distribution', 
                    default='Binomial', 
                    action='store',
                    type='choice',
                    choices=['Binomial', 
                             'Geometric', 
                             'Poisson', 
                             'Binary-0-2',
                             'Binary-0-1-2'],
                    help='distribution to simulate')
  parser.add_option('-o', '--out-dir', 
                    dest='out_dir', 
                    default='result',
                    action='store',
                    type='string',
                    help='output directory')
  parser.add_option('-n', '--nodes', 
                    dest='nodes', 
                    action='callback',
                    type='string',
                    callback=nodes_callback,
                    help='list of the number of nodes')
  parser.add_option('-b', '--basename', 
                    dest='basename', 
                    default='', 
                    action='store',
                    type='string',
                    help='prefix for output files')
  parser.add_option('-p', '',
                    dest='p', 
                    default=0.5, 
                    action='store',
                    type='float',
 help='Probability of success (Binomial,Binary-0-2,Binary-0-1-2,Geometric)')
  parser.add_option('-l', '',
                    dest='l', 
                    default=0.0, 
                    action='store',
                    type='float',
                    help='A new value introduced by Karthik for Binary-0-2')
  parser.add_option('-k', '', 
                    dest='k', 
                    default=2, 
                    action='store',
                    type='int',
                    help='Number of trials (Binomial)')
  parser.add_option('-L', '--lambda', 
                    dest='L', 
                    default=1, 
                    action='store',
                    type='int',
                    help='Mean for Poisson')
  parser.add_option('-t', '--num-experiments', 
                    dest='num_experiments', 
                    default=1, 
                    action='store',
                    type='int',
                    help='Number of experiments to run per tree-size')

  (options, args) = parser.parse_args()

  return options

def main():
  options = parseArguments()

  if not os.path.isdir(options.out_dir):
    print '%s does not exist, creating it ... ' % options.out_dir,
    os.makedirs(options.out_dir)
    print 'done'

  response_file_string = \
'''\
--gen-method %s \
--num-trials %d \
--verbosity 0 \
--print-path true \
--print-tree true \
--measure-mle false \
--test-confidence false \
--lambda %d \
--k %d \
--p %f \
--l %f \
--dump-numbers false \
--use-random-weights false \
--num-threads 1 \
''' % (options.distribution, \
       options.num_experiments, \
       options.L, \
       options.k, \
       options.p,
       options.l)

  if not options.nodes:
    print 'Please enter the number of nodes as a list, e.g., \'-n 10,20,30\''
    return

  for tree_size in options.nodes:

    dist_string = ''
    if options.distribution == 'Binary-0-2':
      dist_string = 'Binary-0-2-p=%f' % options.p
    elif options.distribution == 'Binary-0-1-2':
      dist_string = 'Binary-0-1-2-p=%f' % options.p
    elif options.distribution == 'Binomial':
      dist_string = 'Binomial-p=%f-k=%d' % (options.p, options.k)
    elif options.distribution == 'Geometric':
      dist_string = 'Geometric-p=%f' % options.p
    elif options.distribution == 'Poisson':
      dist_string = 'Poisson-lambda=%d' % options.l

    output_filename = '';
    if options.basename: output_filename = options.basename+'-'+dist_string
    else: output_filename = dist_string

    cur_response_file_string = response_file_string + \
'''\
--n %s \
--dyck-out %s \
--dot-out %s \
''' % (tree_size, \
       tree_size, \
       options.out_dir+'/'+output_filename+'-dyck-n-'+tree_size+'.txt', \
       options.out_dir+'/'+output_filename+'-dot-n-'+tree_size+'.dot')

    print './harness %s' % cur_response_file_string
    os.system('./harness %s' % cur_response_file_string)

# Main starts here
if __name__ == '__main__': main()
