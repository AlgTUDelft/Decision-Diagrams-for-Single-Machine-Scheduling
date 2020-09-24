import itertools
import math
import time
import gc
import heapq

class Solver:
    def __init__(self, dataset_file, limit=3600):
        # Load problem
        self.jobs, self.source, self.sink = extract_jobs(dataset_file)
        
        # Create time windows
        self.time_windows = get_time_windows(sorted(self.jobs))
        self.last_window = self.time_windows[0]
        self.last_partial_job_ids = []
        self.last_prev_job_state = []
        self.priority_queue = []
        self.entry_finder = {}               # mapping of tasks to entries
        self.REMOVED = '<removed-task>'      # placeholder for a removed task
        self.counter = itertools.count()     # unique sequence count

        # Create job state store to avoid duplicates in memory
        self.job_states = {'0': JobState([self.jobs[0].id])}
        
        self.limit=limit
        self.prep_windows()
        
        self.verbose = False

    def slack(self, epsilon):
        return epsilon * self.best_value / (len(self.jobs) - 1)
        
    def solve_fptas(self, epsilon):
        
        self.best_value = 0
        t=time.time()
        timeout = False

        if self.verbose:
            for window in self.time_windows:
                print("{0}".format("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"[window.width]), end="")
            print()

        for key in self.source.path_points:
            for time1, path_point in self.source.path_points[key].items():
                self.add_point(path_point)
        for window in self.time_windows:
            self.y = 0
            self.n = 0
            
            if window.id > 0:
                rem_jobs = [j for j in self.time_windows[window.id - 1].jobs if j not in window.jobs]
                for job in rem_jobs:
                    job.path_points = {}
                    toss = []
                    for key in self.job_states:
                        if job.id in self.job_states[key].job_ids:
                            toss.append(key)
                    for key in toss:
                        del(self.job_states[key])
                if len(rem_jobs) > 0:
                    gc.collect()

            while not self.is_empty():
                if time.time()-t>self.limit:
                    timeout = True
                    break
            
                path_point = self.peek_point()
                job = path_point.job
                job_state = path_point.job_state
                path_point_time = path_point.path_point_time

                if path_point.path_point_time-path_point.job.processing_time >= window.next_start_time:
                    break

                if not self.is_dominated(job, job_state, path_point_time, path_point, window):
                    for successor in job.successors[window.start_time]:
                        if time.time()-t>self.limit:
                             timeout = True
                             break                
                        if successor.id not in job_state.job_ids:
                             new_point,get = self.try_path(job, successor, path_point_time, path_point, job_state, window, self.slack(epsilon))
                             if get:
                                 self.add_point(new_point)

                    #far off
                    for successor in job.far_offs[window.start_time]:
                        if time.time()-t>self.limit:
                             timeout = True
                             break
                        if successor.id not in job_state.job_ids:
                             new_point,get = self.try_path(job, successor, path_point_time, path_point, job_state, window, self.slack(epsilon))
                             if get:
                                self.add_point(new_point)

                else:
                    try:
                        del(job.path_points[job_state.as_key()][path_point_time])
                        del(path_point)
                    except:
                        a=1
                
                self.pop_point()

            if timeout:
                break
            
        return self.sink.get_best_path()
    
    def add_point(self,point):
        if point in self.entry_finder:
            raise Error('point already exists')
        count = next(self.counter)
        entry = [point.point_start_time, count, point]
        self.entry_finder[point] = entry
        heapq.heappush(self.priority_queue, entry)

    def remove_point(self,point):
        try:
            entry = self.entry_finder[point]
            entry[-1] = self.REMOVED
            del self.entry_finder[point]
        except:
            print("remove unexisting point")

    def peek_point(self):
            while self.priority_queue:
                priority, count, point = self.priority_queue[0]
                if point is not self.REMOVED:
                    return point
                else:
                    heapq.heappop(self.priority_queue)
            raise KeyError('peek from an empty priority queue')

    def pop_point(self):
        while self.priority_queue:
            priority, count, point = heapq.heappop(self.priority_queue)
            if point is not self.REMOVED:
                del self.entry_finder[point]
                return point
        raise KeyError('pop from an empty priority queue')

    def is_empty(self):
        while self.priority_queue:
            priority, count, point = self.priority_queue[0]
            if point is not self.REMOVED:
                return 0
            else:
                heapq.heappop(self.priority_queue)
        return 1

    def is_dominated(self, job, job_state, point_time, point, window):
        start_time = point_time
        value = point.value
        if len(job_state.job_ids)>0:
            for subset in job_state.subsets(len(job_state.job_ids)-1):
                subset_key = JobState(subset).as_key()
                if subset_key in job.path_points:
                    for other_point_time, other_point in job.path_points[subset_key].items():
                        if other_point.value >= value and other_point_time <= start_time:
                            return True
        
        if point.previous is not None:
            swapped_job = point.previous.job
            visited = point.must_visits()
            for alt_job in [job for job in window.jobs if job.id not in job_state.job_ids and job not in visited]:
                ids = sorted([j for j in job_state.job_ids if j != swapped_job.id] + [alt_job.id])
                alt_state = JobState(ids)
                alt_state_key = alt_state.as_key()
                if alt_state_key in job.path_points:
                    for other_point_time, other_point in job.path_points[alt_state_key].items():
                        if other_point.value > value and other_point_time <= start_time or other_point.value >= value and other_point_time < start_time:
                            point.must_visit.add(alt_job)
                            break
            for alt_job in [job for job in window.jobs if job.id not in job_state.job_ids and job not in visited and job not in point.must_visit]:
                ids=sorted(job_state.job_ids + [alt_job.id])
                alt_state = JobState(ids)
                alt_state_key = alt_state.as_key()
                if alt_state_key in job.path_points:
                    for other_point_time, other_point in job.path_points[alt_state_key].items():
                        if other_point.value > value and other_point_time <= start_time or other_point.value >= value and other_point_time < start_time:
                            point.must_visit.add(alt_job)
                            break

        must_visits = sorted(list(point.must_visits()), key=lambda j: j.deadline)

        t = start_time
        for i in range(len(must_visits)):
            if t > must_visits[i].latest_start_time:
                self.y += 1

                return True
            t += must_visits[i].processing_time
            t += min([min(mv.setup_times) for mv in must_visits])
        self.n += 1
        return False
    
    """
    Given a job, a successor job, a starting path point and a job state, sees if a new path
    point in the successor job is feasible, not dominated and if so inserts it.
    MUTATING: self.jobs (successor), self.last_window
    """
    def try_path(self, job, successor, path_point_time, path_point, job_state, original_window, slack): 

        temp_point = PathPoint(0, successor, path_point, set(),path_point.VisitedStates,0,path_point.job_state)

        arrival_time = path_point_time
        start_time = max(arrival_time, successor.release_time) + job.setup_time(successor)
        completion_time = start_time + successor.processing_time

        # check feasibility of new path point
        if completion_time > successor.deadline:
            return temp_point,False
        
        # calculate value
        full_value = path_point.value + successor.value
        late_penalty = max(0, successor.penalty_weight * (completion_time - successor.due_date))
        value = full_value - late_penalty
        
        # calculate new job state
        # see if window is the same as last time, saves looking it up
        partial_job_com_ids,partial_job_ids = self.get_partial_job_state(completion_time, successor, job_state, original_window, path_point)
        if successor.latest_start_time>=completion_time:
            partial_job_com_ids.add(successor.id)
            partial_job_ids.append(successor.id)

        job_ids = sorted(partial_job_ids)
 
        new_job_state = JobState(job_ids)
        new_key = new_job_state.as_key()
        
        # Create new state if completely new
        if new_key not in self.job_states:
            self.job_states[new_key] = new_job_state
        
        new_job_state = self.job_states[new_key]

        # create new path point
        new_point =  PathPoint(value, successor, path_point, set(),partial_job_com_ids,completion_time,new_job_state)

        # If the state is new even just in this job, no need to check for domination
        if new_key not in successor.path_points:
            successor.path_points[new_key] = {}
            successor.path_points[new_key][completion_time] = new_point
            return successor.path_points[new_key][completion_time],True

        # check for dominating path_point for the same state (irrespective of windows!)
        for point_time, point in successor.path_points[new_key].items():
            if point_time <= completion_time and point.value + point.max_upper(slack) >= value:
                point.upper = max(point.upper, value - point.value)
                return temp_point,False
        
        # set new best value
        self.best_value = max(self.best_value, value)

        # Delete path points dominated by this new one
        toss = []
        for point_time, point in successor.path_points[new_key].items():
            if point.value <= value + new_point.max_upper(slack) and point_time >= completion_time:
                new_point.upper = max(new_point.upper, point.value - value)
                toss.append(point_time)
        for point_time in toss:
            point = successor.path_points[new_key][point_time]
            self.remove_point(point)
            del(successor.path_points[new_key][point_time])
            del(point)

        successor.path_points[new_key][completion_time] = new_point
        return successor.path_points[new_key][completion_time],True
        
    """
    Finds the time window corresponding to a given time.
    MUTATING: self.last_window
    """
    def get_partial_job_state(self, time, job, job_state, original_window, path_point):
        last_partial_job_com_ids=set()
        self.last_partial_job_ids=[]
        for id in path_point.VisitedStates: 
            if self.jobs[id].latest_start_time>=time:
                last_partial_job_com_ids.add(id)
                if self.jobs[id].latest_start_time>=time+job.setup_time(self.jobs[id]):
                    self.last_partial_job_ids.append(id);
        return last_partial_job_com_ids,self.last_partial_job_ids
    
    """
    Finds potential successors for a job completed in a time window which
    """
    def get_suitable_successors(self, job, window, job_state): #functional
        return [job for job in job.successors[window.start_time] if job.id not in job_state.job_ids]

    def prep_windows(self):
        # Store time window index in each time window
        for i in range(0, len(self.time_windows)):
            self.time_windows[i].id = i
        # Store static info in each job about its time windows and successors
        for window in self.time_windows:
            for job in window.jobs:
                    successors = sorted([successor for successor in self.jobs
                                      # jobs that are reachable
                                      if max(successor.release_time, window.start_time + job.processing_time) + job.setup_time(successor) 
                                      <= successor.latest_start_time
                                      # and are not so far in the future that only the highest value path point is significant
                                      and not min(window.next_start_time - 1, job.latest_start_time) + job.processing_time 
                                      <= successor.release_time
                                     ], key=lambda s: job.setup_time(s))
                    far_offs = sorted([successor for successor in self.jobs
                                    # jobs that are so far in the future that only the highest value path point is significant
                                    if min(window.next_start_time - 1, job.latest_start_time) + job.processing_time 
                                    <= successor.release_time
                                   ], key=lambda s: s.release_time + job.setup_time(s))
                    job.successors[window.start_time] = successors
                    job.far_offs[window.start_time] = far_offs
                    job.windows.append(window)
            self.sink.successors[window.start_time] = []        
    
class Job:
    def __init__(self, i, r, p, e, dd, dl, w, s):
        self.id = i
        self.release_time = r
        self.processing_time = p
        self.value = e
        self.due_date = dd
        self.deadline = dl
        self.penalty_weight = w
        self.setup_times = s
        
        self.latest_start_time = self.deadline - self.processing_time
        
        self.best_path = None
        self.path_points = {} #{job_state.as_key: {start_time: path_point}}
        self.successors = {} #{window.start_time: [Job, sorted by setup time]}
        self.far_offs = {} #same as successors but only those far in the future
        self.windows = [] #list of windows sorted by window start time
        
        self.pretty = ("  [{0}-" + "--{1}-" + "--|{2}|-" + "--{3}]  ({4} - {5}*T)").format(self.release_time, self.processing_time, self.due_date, self.deadline, self.value, self.penalty_weight)

    def __lt__(self, other):
        if self.release_time != other.release_time:
            return self.release_time < other.release_time
        if self.deadline != other.deadline:
            return self.deadline < other.deadline
        return self.id < other.id
    
    def __repr__(self):
        return str(self.id)
    
    def get_job_states(self, size, job_states): #functional
        return [job_states[key] for key in self.path_points if key in job_states and len(job_states[key].job_ids) == size]
    
    def setup_time(self, successor): #functional
        return self.setup_times[successor.id]
    
    def get_path_points(self, job_state, window): #functional
        return [(time, path_point) for time, path_point in self.path_points[job_state.as_key()].items()
                if time-self.processing_time >= window.start_time and time-self.processing_time < window.next_start_time]        
    
    def get_best_path(self):
        value = 0
        path = None
        for job_state in self.path_points:
            for point in self.path_points[job_state].values():
                if point.value > value:
                    path, value = point, point.value
        return path

class Window:
    def __init__(self, start_time, next_start_time, jobs):
        self.start_time = start_time
        self.next_start_time = next_start_time
        self.jobs = jobs
        self.job_ids = sorted([job.id for job in self.jobs])
        self.width = len(jobs)
    
    def __repr__(self):
        return "[{0}\t{1})".format(self.start_time, self.next_start_time)

class JobState:
    __slots__ = ('job_ids', 'key')
    
    def __init__(self, job_ids):
        self.job_ids = job_ids
        self.key = ",".join([str(job_id) for job_id in self.job_ids])
    
    def __len__(self):
        return len(self.job_ids)
    
    def __repr__(self):
        return "(" + self.as_key() + ")"
    
    def as_key(self):
        return self.key
    
    def subsets(self, i):
        return set(itertools.combinations(self.job_ids, i))

class PathPoint:
    __slots__ = ('value', 'job', 'previous', 'must_visit', 'VisitedStates', 'path_point_time','job_state','point_start_time', 'upper')
    
    def __init__(self, value, job, previous, must_visit, VisitedStates, path_point_time, job_state):
        self.value = value
        self.job = job
        self.previous = previous
        self.must_visit = must_visit
        self.VisitedStates = VisitedStates
        self.job_state = job_state
        self.path_point_time = path_point_time
        self.point_start_time= self.path_point_time-job.processing_time
        self.upper = 0

    def max_upper(self, slack_per_job):
        if self.previous is None:
            return slack_per_job - self.upper
        else:
            return self.previous.max_upper(slack_per_job) + slack_per_job - self.upper
    
    def upper_bound(self):
        if self.previous is None:
            return self.upper
        else:
            return self.previous.upper_bound() + self.upper

    def __repr__(self):
        return "(" + str(self.value) + ")"
    
    def get_path(self):
        if self.previous is None:
            return [self.job]
        else:
            return self.previous.get_path() + [self.job]
        
    def must_visits(self):
        if self.previous is None:
            return self.must_visit
        else:
            return self.must_visit | (self.previous.must_visits() - {self.job})

def dataset(n, t, r, i): #functional
    folder = "Dataset_OAS/{0}orders".format(n) + "/Tao{0}".format(t) + "/R{0}".format(r)
    file = "/Dataslack_{0}orders_Tao{1}R{2}_{3}.txt".format(n, t, r, i)
    return folder + file   

def line_to_ints(line): #functional
    return [int(i) for i in line.split(",")]

def line_to_floats(line): #functional
    return [float(i) for i in line.split(",")]

def extract_jobs(dataset): #functional
    file = open(dataset, "r")
    r = line_to_ints(file.readline())
    p = line_to_ints(file.readline())
    dd = line_to_ints(file.readline())
    dl = line_to_ints(file.readline())
    e = line_to_floats(file.readline())
    w = line_to_floats(file.readline())
    s = []

    for i in range(0, len(r)):
        s.append(line_to_ints(file.readline()))

    jobs = []
    for i in range(0, len(r)):
        job = Job(i, r[i], p[i], e[i], dd[i], dl[i], w[i], s[i])
        jobs.append(job)
        
    # Move sink to the end so it won't be inserted inbetween
    jobs[-1].release_time = jobs[-1].deadline
    jobs[-1].due_date = jobs[-1].deadline
    jobs[-1].latest_start_time = jobs[-1].deadline
    jobs[-1].deadline += 1 # algorithm never completes on deadline
        
    # Start path at source
    source_job_state = JobState([jobs[0].id])
    jobs[0].path_points[source_job_state.as_key()] = {}
    source_path_point = PathPoint(jobs[0].value, jobs[0], None, set(), set(),0,source_job_state)
    jobs[0].path_points[source_job_state.as_key()][jobs[0].release_time] = source_path_point
    
    return (jobs, jobs[0], jobs[-1])

def get_time_windows(jobs): #functional
    time_points = []
    for job in jobs:
        time_points.append(job.release_time)
        if job.processing_time>0:
            time_points.append(job.latest_start_time + 1)
    time_points = sorted(list(set(time_points))) # only keep unique entries and sort

    time_windows = []
    for i in range(0, len(time_points) - 1):
        window_start_time = time_points[i]
        next_window_start_time = time_points[i + 1]
        window_jobs = [job for job in jobs if 
                       window_start_time >= job.release_time and 
                       window_start_time <= job.latest_start_time]
        time_window = Window(window_start_time, next_window_start_time, window_jobs)
        time_windows.append(time_window)

    return time_windows


#for phi in [0.05,0.1]:
#    for n in [50,100]:
#        for tau in [9]:
#            for r in [1,3,5,7,9]:
#                for ins in [1,2,3,4,5,6,7,8,9,10]:
#                    if n==100 and r>5:
#                        break
#                    solver = Solver(dataset(n, tau, r, ins),limit=3600)
#                    f = open("oas_FPTAS.txt", 'a+')
#                    t = -time.time()
#                    path = solver.solve_fptas(phi)
#                    t += time.time()

#                    list = []
#                    list.append(path);
#                    point = path.previous
#                    while point!=None:
#                        list.append(point)
#                        point = point.previous

#                    i = len(list)-1;
#                    string1 = "sol:"
#                    while i>=0:
#                        string1 += str(list[i].job.id) +","
#                        i-=1

#                    try:
#                        print("{0},{1},{2},{3}\t{4}\t{5}\t{6}\t{7}".format(n, tau, r, ins, path.value, t,string1,phi), file=f)
#                        print("{0},{1},{2},{3}\t{4}\t{5}\t{6}".format(n, tau, r, ins, path.value, t,phi))
#                    except:
#                        print("{0},{1},{2},{3}\tNo result.".format(n,tau,r,ins)) 

#                    f.close();
#                    del(solver)
#                    del(path)
#                    del(list);
#                    gc.collect()

#for phi in [0.05,0.1]:
#    for v in reversed([v/100 for v in range(0, 102, 2)]):
#        for ins in range(0, 5):
#             filename = "Dataset_OAS/Dataset_bounded_width/CovarianceChange_30orders_c{0}_{1}.txt".format(v, ins)

#             solver = Solver(filename,limit=1800)

#             f = open("cov_FPTAS.txt", 'a+')

#             t = -time.time()
#             path = solver.solve_fptas(phi)
#             t += time.time()

#             try:
#                 list = []
#                 list.append(path);
#                 point = path.previous
#                 while point!=None:
#                     list.append(point)
#                     point = point.previous

#                 i = len(list)-1;
#                 string1 = "sol:"
#                 while i>=0:
#                     string1 += str(list[i].job.id) +","
#                     i-=1

#             #all_t.append(t)
             
#                 print("{0}\tCovarianceChange_{1}_{2}\t{3}\t{4}\t{5}".format(phi,(1-v),ins,path.value, t,string1), file=f)
#                 print("{0}\tCovarianceChange_{1}_{2}\t{3}\t{4}".format(phi,(1-v),ins,path.value,t))
#             except:
#                 print("CovarianceChange_{0}_{1}\tNo results.".format((1-v),ins))

#             del(solver)
#             del(path)
#             del(list);
#             gc.collect()

n = 100

for phi in [0.05,0.1]:
  for R in [1]:
    for w in [3,5,7,9,11,13,15,17,19]:
         for ins in range(0, 5):
             f = open("width_FPTAS.txt", 'a+')
             filename = "Instances_DiffWidth/Dataset_{0}orders_w{1}R{2}_{3}.txt".format(n, w, R, ins)

             solver = Solver(filename,limit=1800)

             t = -time.time()
             path = solver.solve_fptas(phi)
             t += time.time()

             try:
                 print("{0}\t{1},{2},{3}\t{4}\t{5}\t{6}".format(phi,n, w, R, ins, path.value, t))
                 print("{0}\t{1},{2},{3}\t{4}\t{5}\t{6}".format(phi,n, w, R, ins, path.value, t), file=f)
             except:
                 print("{0},{1},{2},{3}\tNo result.".format(n, w, R, ins))  

             del(solver)
             del(path)
             gc.collect()

