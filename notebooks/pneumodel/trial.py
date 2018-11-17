""" This module has functions for working with trial log data produced by simulations containing trial participants. 


"""

from __future__ import division
from bunch import Bunch
from collections import defaultdict
import functools
import json
import numpy as np
import os
import pandas as pd
import scipy.stats
import sys

class Trial:
  """ Contains the simulation results for a trial.

  :param str output_folder: Path to the folder containing output from simulation run(s). 
  :param int run_number: A number between 0 and *n* - 1, where *n* is the number of simulation runs in the output folder.

  """

  def __init__(self, output_folder, run_number):

    self.arms = {}

    # get configuration details
    config_path = os.path.join(output_folder, 'configuration', 'configuration.json')
    with open(config_path) as f:
      config = json.load(f, object_hook=lambda d: Bunch(d))
      trial_config_by_arm = { a.name: a for a in config.simulation.trial_arms}
    
    # get event logs
    path_to_logs = os.path.join(output_folder, 'trial-{}'.format(run_number), 'trial_logs')
    col_names = ['id', 'day', 'event', 'serotype', 'day_of_colonization', 'vaccine_name']
    for f in os.listdir(path_to_logs):
      if f.endswith('.csv'):
        name = f[:-4]
        cfg = trial_config_by_arm[name]
        log = pd.read_csv(os.path.join(path_to_logs, f), header=None, names=col_names, low_memory=False)
        self.arms[name] = TrialArm(cfg, log)

  @property
  def arms(self):
    """ Get the TrialArms of this Trial. """
    return self.__arms

class TrialArm:
  """ Contains the simulation results for one trial arm.

  :param config: A Bunch containing details about the trial arm, with these keys: 

    * name - name of the trial arm
    * num_subjects - number of subjects in the trial arm
    * start_year - the starting year of the trial arm
    * vaccine - the name of the vaccine used in the trial arm
    * schedule - an array of ages (in days) that define the dosing schedule used

  :param log: a Pandas dataframe, where each row describes an event in the trial arm. The data frame contains these columns:

    * id - host id
    * day - day of the event
    * event - type of event (Colonization, Recovery, Vaccination, or Death)
    * serotype - only valid for Colonization and Recovery events
    * day_of_colonization - only valid for *Recovery* events (otherwise, there will be a missing value)
    * vaccine_name - only valid for Vaccination events (otherwise, there will be a missing value)

  .. _Bunch: https://pypi.python.org/pypi/bunch/1.0.1

  """

  def __init__(self, config, log):
    self.config = config
    self.log = log

  @property
  def config(self):
    """ Get the configuration details associated with this TrialArm. """
    return self.__config

  @property
  def log(self):
    """ Get the event log associated with this TrialArm. """
    return self.__log

def get_num_runs(folder):
  n = 0
  while(True):
    path = os.path.join(folder, 'trial-{}'.format(n))
    if not os.path.exists(path):
      break
    else:
      n += 1
  return n

def load_trials(folder):
  return [Trial(folder, i) for i in range(get_num_runs(folder))]


def get_incidence(arm, fup_days, serotype=None, first_day=None):
  """ Find the incidence of colonization in a trial arm over a follow-up period. (Currently assumes no deaths over follow-up period.)

  :param arm: A :class:`pneumodel.trial.TrialArm` instance
  :param fup_days: The length of the follow-up period, in days.
  :param serotype: The serotypes of interest. If None, then all serotypes are included.
  :param first_day: The first day of the follow-up period, relative to the start of the trial (i.e. the birth of the trial arm subjects). If None, then the first day is the day the last vaccine dose is administered.

  :return float: The incidence rate (number of events / person-days) 

  .. _Bunch: https://pypi.python.org/pypi/bunch/1.0.1
  """

  # set the first day of follow up
  if first_day is None:
    first_day = arm.config.start_year * 365 + max(arm.config.schedule)
  else:
    first_day = arm.config.start_year * 365 + first_day
  
  # an incident is a colonization event within the followup period
  is_incident = (arm.log.event == 'C') & (arm.log.day >= first_day) & (arm.log.day <= first_day + fup_days)
  if serotype:
    is_incident = is_incident & (arm.log.serotype == serotype)

  # person-day calculation
  person_days = arm.config.num_subjects * fup_days # assumes no one dies
  is_death = (arm.log.event == 'D') & (arm.log.day >= first_day) & (arm.log.day <= first_day + fup_days)
  deaths = arm.log[is_death]
  lost_days = np.sum((first_day + fup_days) - deaths.day)
  person_days -= lost_days

  # incidence rate is number of incidents / person-days
  num_incidents = np.sum(is_incident)
  return num_incidents / person_days

def get_mean_duration(arm, fup_days, serotype=None, first_day=None):
  """ Find the mean duration of colonizations in a trial arm over a follow-up period.

  :param arm: A :class:`pneumodel.trial.TrialArm` instance
  :param fup_days: The length of the follow-up period, in days.
  :param serotype: The serotypes of interest. Specify None to include all serotypes.
  :param first_day: The first day of the follow-up period, relative to the start of the trial (i.e. the birth of the trial arm subjects). If None, then the first day is the day the last vaccine dose is administered.

  .. _Bunch: https://pypi.python.org/pypi/bunch/1.0.1

  """

  # set the first day of follow up
  if first_day is None:
    first_day = arm.config.start_year * 365 + max(arm.config.schedule)
  else:
    first_day = arm.config.start_year * 365 + first_day
  
  # get all the recovery events within our follow-up period, whose day of colonization occurred after vaccination
  is_recovery = (arm.log.event == 'R') & (arm.log.day_of_colonization >= first_day) & (arm.log.day < first_day + fup_days)
  if serotype:
    is_recovery = is_recovery & (arm.log.serotype == serotype)
  
  # calculate durations from information logged in the recovery events
  recoveries = arm.log[is_recovery]
  durations = recoveries.day - recoveries.day_of_colonization
  return np.mean(durations)


def _get_num_colonized(hosts, serotypes=None):
  """ Calculates the frequency of colonization.

  :param hosts: A dictionary mapping host IDs to a list of serotypes carried by that host.
  :param serotypes : Serotypes of interest. If None, all serotypes are considered.

  """
  count = 0
  for host_id, carriage in hosts.iteritems():
    if serotypes:
      is_colonized = any(s in carriage for s in serotypes)
    else:
      is_colonized = len(carriage) > 0
    
    if is_colonized:
      count += 1

  return count

def get_prevalence(arm, sampling_days, serotypes=None, first_day=None):
  """ A generator that yields prevalence of colonization in a trial arm on different sampling days.

  :param arm: A :class:`pneumodel.trial.TrialArm` instance
  :param sampling_days: Days to sample the prevalence, relative to the first day of the follow up period.
  :param serotypes: The serotypes of interest. Specify None to include all serotypes.
  :param first_day: The first day of the follow-up period, relative to the start of the trial (i.e. the birth of the trial arm subjects). If None, then the first day is the day the last vaccine dose is administered.

  :yield: The point prevalence for each sampling day.

  .. _Bunch: https://pypi.python.org/pypi/bunch/1.0.1

  """
  
  # set the first day of follow up
  if first_day is None:
    first_day = arm.config.start_year * 365 + max(arm.config.schedule)
  else:
    first_day = arm.config.start_year * 365 + first_day

  sampling_days = sorted(sampling_days)
  sampling_days = (first_day + d for d in sampling_days)
  
  # for finding our place temporally
  previous_day = 0
  next_sampling_day = next(sampling_days)
  
  # each "host" is a list of serotypes
  hosts = defaultdict(list)

  # go through all the events
  for index, row in arm.log.iterrows():

    # if the event is the first one after a sampling day, report the prevalence we have so far
    if row.day >= next_sampling_day and previous_day < next_sampling_day:
      yield float(_get_num_colonized(hosts, serotypes)) / arm.config.num_subjects # assumes no ones dies during trial
      try:
        next_sampling_day = next(sampling_days)
      except StopIteration:
        return
    
    # otherwise, process the event
    try:
      if row.event == 'C':
        hosts[row.id].append(row.serotype)
      elif row.event == 'R':
        hosts[row.id].remove(row.serotype)
      elif row.event == 'D':
        if hosts.has_key(row.id):
          del hosts[row.id]
    except:
      print arm.config
      print row
      print hosts[row.id]
      raise RuntimeError('oops')
    
    previous_day = row.day

def get_prevalence_trajectory(arm, sampling_days, serotypes=None, first_day=None):
  """ Returns a list of colonization prevalences in a trial arm on different sampling days.

  :param sampling_days: Days to sample the prevalence, relative to the first day of the follow up period.
  :param serotypes: The serotypes of interest. Specify None to include all serotypes.
  :param first_day: The first day of the follow-up period, relative to the start of the trial (i.e. the birth of the trial arm subjects). If None, then the first day is the day the last vaccine dose is administered.

  :return: A list of point prevalences, one for each sampling day.

  .. _Bunch: https://pypi.python.org/pypi/bunch/1.0.1

  """
  return list(get_prevalence(arm, sampling_days, serotypes, first_day))


def sample_size(p1, p2, a=0.05, b=0.8):
  """ Calculates the sample size needed for a two-sample test of equal proportions (see link_). Assumes equal sample sizes.

  :param p1: Proportion in one sample
  :param p2: Proportion in the other sample
  :param a: Type I error probability
  :param b: Desired power

  :returns: The total sample size needed (size of sample 1 + size of sample 2).

  .. _link: http://hansheng.gsm.pku.edu.cn/pdf/2007/prop.pdf
  """
  p1 = np.array(p1)
  p2 = np.array(p2)
  z_half_a = scipy.stats.norm.ppf(1 - a / 2.)
  z_b      = scipy.stats.norm.ppf(b)
  e        = p1 - p2
  return 2 * (((z_half_a + z_b)**2) / (e**2)) * (p1 * (1 - p1) + p2 * (1 - p2))


def evaluate_across_runs(func, trials, arm_name, *args, **kwargs):
  """ Executes a function on the same trial arm over multiple simulation runs.

  :param func: This function should take as its first argument a :class:`pneumodel.trial.TrialArm`
  :param trials: An iterable of :class:`pneumodel.trial.Trial` instances 
  :param arm_name: Name of the trial arm of interest.
  :param args: Any positional arguments for func.
  :param kwargs: Any keyword arguments for func.

  :returns: A numpy array containing the results of func evaluated on each simulation run.

  """
  return np.array([func(trial.arms[arm_name], *args, **kwargs) for trial in trials])


##########################
### UPDATES 03/06/2017 ###
##########################

def load_trajectories_single_trial(run_path, index_by):
  trj = {}
  
  # navigate to the folder with the relevant data
  path_to_logs = os.path.join(run_path, 'trial_prevalence')

  # check each file
  for f in os.listdir(path_to_logs):
    path = os.path.join(path_to_logs, f)
    
    # if it's a csv file (all files in there should be)
    if os.path.isfile(path) and f.endswith('.csv'):
      
      # read in data and calculate prevalence
      arm_name = f.split('.csv')[0]
      df = pd.read_csv(path, header=None)
      df.columns = ['dsb', 'dsv', 'num_alive', 'num_colonized']
      df['prevalence'] = df.num_colonized / df.num_alive
      df.drop('num_alive', axis=1, inplace=True)
      df.drop('num_colonized', axis=1, inplace=True)
      
      # select calendar to index by
      if index_by == 'birth':
        df.drop('dsv', axis=1, inplace=True)
        df.set_index('dsb', inplace=True)
      elif index_by == 'last_dose':
        df.drop('dsb', axis=1, inplace=True)
        df.set_index('dsv', inplace=True)
      else:
        raise ValueError('index_by should be "birth" or "last_dose"')
        
      # add this dataframe to our dictionary
      trj[arm_name] = df
      
  return trj

def load_trajectories(output_path, index_by, max_runs=None):
  # determine number of runs
  num_trials = get_num_runs(output_path)
  assert(num_trials) > 0

  # load trajectories for every run
  paths = [os.path.join(output_path, 'trial-{}'.format(i)) for i in range(num_trials)]
  trajectories = [load_trajectories_single_trial(p, index_by) for p in paths]
  
  # collect results across runs by trial arm
  result = {}
  for arm_name in trajectories[0].keys():
    # merge dataframes from each run
    dfs = [trj[arm_name] for trj in trajectories]
    merged = pd.concat(dfs, axis=1)
    
    # rename the columns
    merged.columns = ['run-{}'.format(i) for i in range(num_trials)]
    result[arm_name] = merged.iloc[:,:max_runs]

  return result

