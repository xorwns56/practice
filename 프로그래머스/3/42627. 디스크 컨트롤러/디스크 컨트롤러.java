import java.util.*;
class Job{
    int number;
    int time;
    int duration;
    Job(int number, int time, int duration){
        this.number = number;
        this.time = time;
        this.duration = duration;
    }
}
class Solution {
    public int solution(int[][] jobs) {
        List<Job> job_list = new ArrayList<>();
        for(int i = 0; i < jobs.length; i++) job_list.add(new Job(i, jobs[i][0], jobs[i][1]));
        Collections.sort(job_list, (job1, job2)->{ return job1.time - job2.time; });   
        PriorityQueue<Job> q = new PriorityQueue<Job>((job1, job2)->{
            if(job1.duration != job2.duration) return job1.duration - job2.duration;
            if(job1.time != job2.time) return job1.time - job2.time;
            return job1.number - job2.number;
        });
        Job curr = null;
        int time = 0;
        int count = 0;
        int sum = 0;
        while(count < jobs.length){
            while(job_list.size() > 0 && time >= job_list.get(0).time) q.add(job_list.remove(0));
            if(curr == null && q.size() > 0){
                curr = q.remove();
                time += curr.duration;
                sum += time - curr.time;
                curr = null;
                count++;
            }else time++;
        }
        return sum / count;
    }
}