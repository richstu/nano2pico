#!/bin/env python3
import os
import subprocess
import sys

DEBUG = False

# envDict: { key: value }
def findEnviornment(scriptname, envDict):
  if not os.path.isfile(scriptname):
    print ("[Error] Can't find script:"+scriptname)

  command = ['env', '-i', 'bash', '-c', 'source '+scriptname+' && env']
  proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell=True)
  for line in proc.stdout:
    if '{' in line.decode('utf-8'): continue
    if '}' in line.decode('utf-8'): continue
    t_array= line.decode('utf-8').split("=",1)
    key=t_array[0]
    value=t_array[1].rstrip('\n')
    envDict[key] = value
  proc.communicate()

def returnEnviornment(scriptname):
  envDict = {}
  findEnviornment(scriptname, envDict)
  return envDict

def addRootEnv(_env):
  _env.Append (CCFLAGS = '-isystem `root-config --incdir`' )
  if DEBUG:
    _env.Append (CCFLAGS = '-g' )
  _env.Append (CCFLAGS = '`root-config --cflags`' )
  _env.Append (LINKFLAGS = '`root-config --glibs`') 
  _env.Append (LINKFLAGS = '`root-config --ldflags`')
  _env.Append (LINKFLAGS = ['-lRooFit', '-lGenVector', '-lRooStats', '-lRooFitCore'])

def addWarningEnv(_env):
  _env.Append (CCFLAGS = ['-pedantic', 
                          '-Wall', '-Wextra', '-Werror', '-Wold-style-cast', 
                          '-Wcast-align', '-Wcast-qual', '-Wdisabled-optimization', 
                          '-Wformat=2', '-Wformat-nonliteral', '-Wformat-security', 
                          '-Wformat-y2k', '-Winit-self', '-Winvalid-pch', '-Wlong-long', 
                          '-Wmissing-format-attribute', '-Wmissing-include-dirs',
                          '-Wpacked', '-Wpointer-arith', '-Wredundant-decls', '-Wstack-protector', 
                          '-Wundef', '-Wvariadic-macros', '-Wmissing-noreturn', 
                          '-Wwrite-strings', '-Wctor-dtor-privacy', '-Wnon-virtual-dtor', '-Wsign-promo', '-Wsign-compare', 
                          '-Wunreachable-code', 
                          '-Woverloaded-virtual', '-Wshadow', '-Wswitch-default', '-Wswitch-enum', '-Wunused', 
                          #'-Wsign-conversion', '-Wfloat-equal', '-Wunsafe-loop-optimizations', 
                         ])

def addExternalEnv(_env):
  _env.Append (CCFLAGS = '-isystem external_inc -std=c++17' )

def addBasicEnv(_env):
  if not DEBUG:
    _env.Append (CCFLAGS = '-O2')

def addKernelEnv(_env):
  _env['kernel'] = getKernel()

def getKernel():
  return subprocess.check_output("uname -r | cut -d '-' -f1", shell=True, universal_newlines=True).rstrip()


SConsignFile('kernel/'+getKernel()+'/sconsign.dblite')

# [Requried] Export SET_ENV_PATH to script that sets root and scons environment
analysisEnv = Environment(ENV = returnEnviornment(os.environ['SET_ENV_PATH']))

addBasicEnv(analysisEnv)
addKernelEnv(analysisEnv)
addRootEnv(analysisEnv)
addExternalEnv(analysisEnv)
addWarningEnv(analysisEnv)

exportEnv = analysisEnv
SConscript('SConscript', variant_dir='build/'+analysisEnv['kernel'], duplicate=0, exports="exportEnv")
