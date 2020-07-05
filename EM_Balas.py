import itertools
import math
import time
# import numpy
import gc
#import xlsxwriter


class Solver:
    def __init__(self, dataset_file, w=100, heuristic=False, heuristic_function=None):
        # Load problem
        self.jobs, self.source, self.sink = extract_jobs(dataset_file,heuristic_function)

        # Create time windows
        self.time_windows = get_time_windows(sorted(self.jobs))
        self.last_window = self.time_windows[0]
        self.last_partial_job_ids = []
        self.last_prev_job_state = []

        # Create job state store to avoid duplicates in memory
        self.job_states = {'0,1': JobState([self.jobs[0].id], 1)}

        self.heuristic = heuristic

        self.w = w
        self.w2 = w

        self.heuristic_function = heuristic_function

        self.prep_windows()

        self.verbose = False

    def solve_exact(self):
        t = time.time()
        timeout = False
        if self.verbose:
            for window in self.time_windows:
                print("{0}".format("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"[window.width]),
                      end="")
            print()
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
                        del (self.job_states[key])
                if len(rem_jobs) > 0:
                    gc.collect()

            for job_state_size in range(0, window.width + 1):
                for job in window.jobs:
                    job_states = job.get_job_states(job_state_size, self.job_states)
                    for job_state in job_states:
                        path_points = job.get_path_points(job_state, window)
                        if len(path_points) > 0:
                            max_value = max([other.value for other_time, other in path_points])
                            max_path_time, max_path = \
                                [(point_time, point) for point_time, point in path_points if point.value == max_value][
                                    0]
                        for path_point_time, path_point in path_points:

                            if time.time() - t > 1800:  # or psutil.virtual_memory().percent > 90:
                                timeout = True
                                break
                            if not self.is_dominated(job, job_state, path_point_time, path_point, window):
                                for successor in job.successors[window.start_time]:
                                    if time.time() - t > 1800:  # or psutil.virtual_memory().percent > 90:
                                        timeout = True
                                        break
                                    if successor.id not in job_state.job_ids:
                                        self.try_path(job, successor, path_point_time, path_point, job_state, window)

                                # Only consider last path point for far off vertices
                                if path_point.value == max_value:
                                    for successor in job.far_offs[window.start_time]:
                                        if time.time() - t > 1800:  # or psutil.virtual_memory().percent > 90:
                                            timeout = True
                                            break
                                        if successor.id not in job_state.job_ids:
                                            self.try_path(job, successor, max_path_time, max_path, job_state, window)
                            else:
                                del (job.path_points[job_state.as_key()][path_point_time])
                                del (path_point)
                            if timeout:
                                break
                        if timeout:
                            break
                    if timeout:
                        break
                if timeout:
                    break
            if timeout:
                break

            if self.verbose:
                if self.y >= self.n:
                    print("y", end="")
                else:
                    print(min(9, int(self.n / (self.y + 1))), end="")
        # memoryout = False
        # if psutil.virtual_memory().percent > 90:
        #    memoryout=True
        return self.sink.get_best_path()

    def is_dominated(self, job, job_state, point_time, point, window):
        start_time = point_time
        value = point.value
        if len(job_state.job_ids) > 0:
            for subset in job_state.subsets(len(job_state.job_ids) - 1):
                subset_key = JobState(subset, len(subset)).as_key()
                if subset_key in job.path_points:
                    for other_point_time, other_point in job.path_points[subset_key].items():
                        if other_point.value >= value and other_point_time <= start_time:  # and subset_key!= job_state.key:
                            return True

        if point.previous is not None:
            swapped_job = point.previous.job
            visited = point.must_visits()
            for alt_job in [job for job in window.jobs if job.id not in job_state.job_ids and job not in visited]:
                ids = sorted([j for j in job_state.job_ids if j != swapped_job.id] + [alt_job.id])
                alt_state = JobState(ids, len(ids))
                alt_state_key = alt_state.as_key()
                if alt_state_key in job.path_points:
                    for other_point_time, other_point in job.path_points[alt_state_key].items():
                        if other_point.value > value and other_point_time <= start_time or other_point.value >= value and other_point_time < start_time:
                            point.must_visit.add(alt_job)
                            break
            for alt_job in [job for job in window.jobs if
                            job.id not in job_state.job_ids and job not in visited and job not in point.must_visit]:
                ids = sorted(job_state.job_ids + [alt_job.id])
                alt_state = JobState(ids, len(ids))
                alt_state_key = alt_state.as_key()
                if alt_state_key in job.path_points:
                    for other_point_time, other_point in job.path_points[alt_state_key].items():
                        if other_point.value > value and other_point_time <= start_time or other_point.value >= value and other_point_time < start_time:
                            point.must_visit.add(alt_job)
                            break

        must_visits = sorted(list(point.must_visits()), key=lambda j: j.deadline)

        t = start_time  # + job.processing_time
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
    TODO: optimise, most used method
    MUTATING: self.jobs (successor), self.last_window
    """

    def try_path(self, job, successor, path_point_time, path_point, job_state, original_window):

        #apply balas neighbourhood
        if successor.index + self.w <= path_point.index: # successor must be scheduled before path_point
            return 

        arrival_time = path_point_time  # + job.processing_time
        start_time = max(arrival_time, successor.release_time) + job.setup_time(successor)
        completion_time = start_time + successor.processing_time

        # check feasibility of new path point
        if completion_time > successor.deadline:
            return

            # calculate value
        full_value = path_point.value + successor.value
        late_penalty = max(0, successor.penalty_weight * (completion_time - successor.due_date))
        value = full_value - late_penalty

        # calculate new job state
        # see if window is the same as last time, saves looking it up
        partial_job_com_ids, partial_job_ids = self.get_partial_job_state(completion_time, successor, job_state,
                                                                          original_window, path_point)
        if successor.latest_start_time >= completion_time:
            partial_job_com_ids.add(successor.id)
            partial_job_ids.append(successor.id)

        job_ids = sorted(partial_job_ids)

        size = 0

        if start_time < original_window.next_start_time and len(job_ids) <= job_state.size:
            size = job_state.size + 1
        else:
            size = len(job_ids)

        new_job_state = JobState(job_ids, size)
        new_key = new_job_state.as_key()
        # new_job_state.job_com_ids = partial_job_com_ids

        #Calculate current max index
        index = max(path_point.index,successor.index)

        # Create new state if completely new
        if new_key not in self.job_states:
            self.job_states[new_key] = new_job_state
        # else:
        #    self.job_states[new_key].job_com_ids= list(set(new_job_state.job_com_ids).union(set(self.job_states[new_key].job_com_ids)))

        new_job_state = self.job_states[new_key]

        # If the state is new even just in this job, no need to check for domination
        if new_key not in successor.path_points:
            successor.path_points[new_key] = {}
            successor.path_points[new_key][completion_time] = PathPoint(value, successor, path_point, set(),
                                                                        partial_job_com_ids,index)
            return

        # check for dominating path_point for the same state (irrespective of windows!)
        for point_time, point in successor.path_points[new_key].items():
            # if point_time == completion_time:
            #    self.job_states[new_key].job_com_ids= partial_job_com_ids|self.job_states[new_key].job_com_ids
            if point_time <= completion_time and point.value >= value:
                return

        # Delete path points dominated by this new one
        toss = []
        for point_time, point in successor.path_points[new_key].items():
            if point.value <= value and point_time >= completion_time:
                toss.append(point_time)
        for point_time in toss:
            point = successor.path_points[new_key][point_time]
            del (successor.path_points[new_key][point_time])
            del (point)
        # for point_time, point in successor.path_points[new_key].items():

        # if value<successor.path_points[new_key][completion_time].value:
        #    aaaa=1;

        successor.path_points[new_key][completion_time] = PathPoint(value, successor, path_point, set(),
                                                                    partial_job_com_ids,index)
        return

    """
    Finds the time window corresponding to a given time.
    TODO: optimise, most used method
    MUTATING: self.last_window
    """

    def get_partial_job_state(self, time, job, job_state, original_window, path_point):
        # track last window because it is not all that unlikely that this remains the same for many path points
        # if time < self.last_window.start_time or time >= self.last_window.next_start_time or job_state.job_ids != self.last_prev_job_state:
        #    if time >= self.last_window.next_start_time:
        #        for i in range(self.last_window.id + 1, len(self.time_windows)):
        #            if time < self.time_windows[i].next_start_time:
        #                self.last_window = self.time_windows[i]
        #                break
        #    else:
        #        for i in range(original_window.id, self.last_window.id):
        #            if time < self.time_windows[i].next_start_time:
        #                self.last_window = self.time_windows[i]
        #                break
        #    self.last_prev_job_state = job_state.job_ids
        #    self.last_partial_job_ids = [i for i in job_state.job_ids if i in self.last_window.job_ids]
        last_partial_job_com_ids = set()
        self.last_partial_job_ids = []
        # self.jobs
        # for Job in original_window.jobs:
        #    if Job.id in job_state.job_ids and Job.latest_start_time>=time: path_point.VisitedStates: 
        for id in path_point.VisitedStates:
            if self.jobs[id].latest_start_time >= time:
                last_partial_job_com_ids.add(id)
                if self.jobs[id].latest_start_time >= time + job.setup_time(self.jobs[id]) and abs(self.jobs[id].index -job.index) < self.w:
                    self.last_partial_job_ids.append(id);
        return last_partial_job_com_ids, self.last_partial_job_ids

    """
    Finds potential successors for a job completed in a time window which
    TODO: optimise, most used method
    """

    def get_suitable_successors(self, job, window, job_state):  # functional
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
                                     # F(i, t_i)
                                     if max(successor.release_time,
                                            window.start_time + job.processing_time) + job.setup_time(successor)
                                     <= successor.latest_start_time
                                     # and are not so far in the future that only the highest value path point is significant
                                     and not min(window.next_start_time - 1,
                                                 job.latest_start_time) + job.processing_time
                                             <= successor.release_time
                                             # Apply balas neighborhood 
                                     and (abs(successor.index -job.index) < self.w or successor.index==len(self.jobs)-1)
                                     ], key=lambda s: job.setup_time(s)
                                                                       )
                far_offs = sorted([successor for successor in self.jobs
                                   # jobs that are so far in the future that only the highest value path point is significant
                                   if min(window.next_start_time - 1, job.latest_start_time) + job.processing_time
                                   <= successor.release_time
                                   # Apply balas neighborhood 
                                     and (abs(successor.index -job.index) < self.w or successor.index==len(self.jobs)-1)
                                   ], key=lambda s: s.release_time + job.setup_time(s)
                                    )
                #if self.heuristic:
                #    successors = self.apply_heuristic(successors, self.w, self.heuristic_function)
                #    far_offs = self.apply_heuristic(far_offs, self.w2, self.heuristic_function)

                job.successors[window.start_time] = successors
                job.far_offs[window.start_time] = far_offs
                job.windows.append(window)
            self.sink.successors[window.start_time] = []

    def apply_heuristic(self, job_array, heuristic_size, heuristic_function):
        if len(job_array) < heuristic_size:
            return job_array
        new_array = sorted(job_array, key=heuristic_function)[:heuristic_size + 1]
        return new_array


def latest_start_time(x):
    return x.latest_start_time


def value(x):
    return -x.value


def unit_value(x):
    if x.processing_time != 0:
        return -x.value / x.processing_time
    return 0


def processing_time(x):
    return x.processing_time


def combined(x):
    if x.value != 0:
        return x.latest_start_time ** 2 / x.value
    return x.latest_start_time ** 2


class Job:
    # why not also use slots here?
    def __init__(self, i, r, p, e, dd, dl, w, s):
        self.id = i
        self.release_time = r
        self.processing_time = p
        self.value = e
        self.due_date = dd
        # w * (dl - dd) <= e
        self.deadline = dl
        self.penalty_weight = w
        # array
        self.setup_times = s

        self.index = None

        self.latest_start_time = self.deadline - self.processing_time

        self.best_path = None
        self.path_points = {}  # {job_state.as_key: {start_time: path_point}}
        self.successors = {}  # {window.start_time: [Job, sorted by setup time]}
        self.far_offs = {}  # same as successors but only those far in the future
        self.windows = []  # list of windows sorted by window start time

        self.pretty = ("  [{0}-" + "--{1}-" + "--|{2}|-" + "--{3}]  ({4} - {5}*T)").format(self.release_time,
                                                                                           self.processing_time,
                                                                                           self.due_date, self.deadline,
                                                                                           self.value,
                                                                                           self.penalty_weight)

    def __lt__(self, other):
        if self.release_time != other.release_time:
            return self.release_time < other.release_time
        if self.deadline != other.deadline:
            return self.deadline < other.deadline
        return self.id < other.id

    def __repr__(self):
        return str(self.id)

    def get_job_states(self, size, job_states):  # functional
        # for key in self.path_points:
        #    if key in job_states:
        #        if job_states[key].size==size:
        #            aaaa=1
        return [job_states[key] for key in self.path_points if key in job_states and job_states[key].size == size]

    def setup_time(self, successor):  # functional
        return self.setup_times[successor.id]

    def get_path_points(self, job_state, window):  # functional
        return [(time, path_point) for time, path_point in self.path_points[job_state.as_key()].items()
                if
                time - self.processing_time >= window.start_time and time - self.processing_time < window.next_start_time]

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
    __slots__ = ('job_ids', 'key', 'size')

    def __init__(self, job_ids, size):
        self.job_ids = job_ids
        self.key = ",".join([str(job_id) for job_id in self.job_ids])
        self.key += "," + str(size)
        self.size = size
        # self.job_com_ids=set(); ,'job_com_ids'

    def __len__(self):
        return self.size

    def __repr__(self):
        return "(" + self.as_key() + ")"

    def as_key(self):
        return self.key

    def subsets(self, i):
        return set(itertools.combinations(self.job_ids, i))


class PathPoint:
    __slots__ = ('value', 'job', 'previous', 'must_visit', 'VisitedStates', 'index')

    def __init__(self, value, job, previous, must_visit, VisitedStates, index):
        self.value = value
        self.job = job
        self.previous = previous
        self.must_visit = must_visit
        self.VisitedStates = VisitedStates
        self.index = index # Current max index 

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


def dataset(n, t, r, i):  # functional
    folder = "Dataset_OAS/{0}orders".format(
        n) + "/Tao{0}".format(t) + "/R{0}".format(r)
    file = "/Dataslack_{0}orders_Tao{1}R{2}_{3}.txt".format(n, t, r, i)
    return folder + file


def line_to_ints(line):  # functional
    return [int(i) for i in line.split(",")]


def line_to_floats(line):  # functional
    return [float(i) for i in line.split(",")]


def extract_jobs(dataset,heuristic_function):  # functional
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
    jobs[-1].deadline += 1  # algorithm never completes on deadline

    # Calculate index
    index_jobs = [];
    index_jobs = sorted(jobs, key=heuristic_function)
    for i in range(0,len(index_jobs)):
        index_jobs[i].index = i

    # Start path at source
    source_job_state = JobState(job_ids=[jobs[0].id], size=1)

    jobs[0].path_points[source_job_state.as_key()] = {}

    source_path_point = PathPoint(value=jobs[0].value,
                                  job=jobs[0],
                                  previous=None,
                                  must_visit=set(),
                                  VisitedStates=set(),
                                  index=0)
    jobs[0].path_points[source_job_state.as_key()][jobs[0].release_time] = source_path_point

    return (jobs, jobs[0], jobs[-1])


def get_time_windows(jobs):  # functional
    time_points = []
    for job in jobs:
        time_points.append(job.release_time)
        if job.processing_time > 0:
            time_points.append(job.latest_start_time + 1)
    time_points = sorted(list(set(time_points)))  # only keep unique entries and sort

    time_windows = []
    for i in range(0, len(time_points) - 1):
        window_start_time = time_points[i]
        next_window_start_time = time_points[i + 1]
        window_jobs = [job for job in jobs if
                       job.release_time <= window_start_time <= job.latest_start_time]
        time_window = Window(window_start_time, next_window_start_time, window_jobs)
        time_windows.append(time_window)

    return time_windows


def filename(n, tau, r, i):
    return "Dataset_OAS/{0}orders/Tao{1}/R{2}/Dataslack_{0}orders_Tao{1}R{2}_{3}.txt".format(n, tau, r, i)

def width_filename(width, r, ins):
    return "Dataset_OAS/Changing width/Dataset_100orders_w{0}R{1}_{2}.txt".format(width, r, ins)



#heuristic_function = latest_start_time


#for n in [50,100]:
#        for tau in [9]:
#            for r in [1,3,5,7,9]:
#                for ins in [1,2,3,4,5,6,7,8,9,10]:
#                  for w in [5,7,9,11,13,15,17,19]:
#                    solver = Solver(dataset(n, tau, r, ins), heuristic=True, w=w, heuristic_function=heuristic_function)
#                    f = open("resultwithbalas.txt", 'a+')
#                    t = -time.time()
#                    path = solver.solve_exact()
#                    t += time.time()

#                    try:
#                        list = []
#                        list.append(path);

#                        point = path.previous
#                        while point!=None:
#                            list.append(point)
#                            point = point.previous

#                        i = len(list)-1;
#                        string1 = "sol:"
#                        while i>=0:
#                            string1 += str(list[i].job.id) +","
#                            i-=1

                    
#                        print("{0}\t{1},{2},{3}\t{4}\t{5}\t{6}\t{7}".format(w,n, tau, r, ins, path.value, t,string1), file=f)
#                        print("{0}\t{1},{2},{3}\t{4}\t{5}\t{6}".format(w,n, tau, r, ins, path.value, t))
#                    except:
#                        print("{0},{1},{2},{3}\tNo result.".format(n,tau,r,ins), file=f) 
                    

#                    f.close();
#                    del(solver)
#                    del(path)
#                    del(list);
#                    gc.collect()
                
#heuristic_function = latest_start_time     
#for wh in [7,9,11,13,15,17,19]:
#   for v in reversed([v/100 for v in range(0, 102, 2)]):
#        for ins in range(0, 5):
#             filename = "Dataset_OAS/Changing covariance/CovarianceChange_30orders_c{0}_{1}.txt".format(v, ins)

#             solver = Solver(filename, heuristic=True, w=wh, heuristic_function=heuristic_function)

#             f = open("resultwithbalas.txt", 'a+')

#             t = -time.time()
#             path = solver.solve_exact()
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
             
#                 print("{0}\tCovarianceChange_{1}_{2}\t{3}\t{4}\t{5}".format(wh,(1-v),ins,path.value, t,string1), file=f)
#                 print("{0}\tCovarianceChange_{1}_{2}\t{3}\t{4}".format(wh,(1-v),ins,path.value,t))
#             except:
#                 print("CovarianceChange_{0}_{1}\tNo results.".format((1-v),ins), file=f)

#             del(solver)
#             del(path)
#             del(list);
#             gc.collect()

heuristic_function = latest_start_time
n = 100

for wh in [5,7,9,11,13,15,17,19]:
  for R in [9]:
    for w in [3,5,7,9,11,15,17,19]:
         for ins in range(0, 10):
             f = open("resultwithbalas.txt", 'a+')
             filename = "Dataset_OAS/Changing width/Dataset_{0}orders_w{1}R{2}_{3}.txt".format(n, w, R, ins)

             solver = Solver(filename, heuristic=True, w=wh, heuristic_function=heuristic_function)

             

             t = -time.time()
             path = solver.solve_exact()
             t += time.time()

             #all_t.append(t)
             try:
                 print("{0}\t{1},{2},{3}\t{4}\t{5}\t{6}".format(wh,n, w, R, ins, path.value, t))
                 print("{0}\t{1},{2},{3}\t{4}\t{5}\t{6}".format(wh,n, w, R, ins, path.value, t), file=f)
             except:
                 print("{0},{1},{2},{3}\tNo result.".format(n, w, R, ins), file=f)  

             del(solver)
             del(path)
             gc.collect()