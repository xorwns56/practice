import java.util.*;
class Stage implements Comparable<Stage>{
    int index;
    int user;
    int fail_user;
    Stage(int index, int user, int fail_user){
        this.index = index;
        this.user = user;
        this.fail_user = user;
    }
    @Override
    public int compareTo(Stage s){
        double fail_rate = this.user != 0 ? (double)this.fail_user / this.user : 0;
        double fail_rate2 = s.user != 0 ? (double)s.fail_user / s.user : 0;
        if(fail_rate == fail_rate2) return this.index - s.index;
        return fail_rate > fail_rate2 ? -1 : 1;
    }
}
class Solution {
    public int[] solution(int N, int[] stages) {
        Stage[] stage = new Stage[N];
        for(int i = 0; i < stages.length; i++){
            for(int j = 0; j < N; j++){
                if(i == 0) stage[j] = new Stage(j + 1, 0, 0);
                if(j + 1 <= stages[i]){
                    stage[j].user++;
                    if(j + 1 == stages[i]) stage[j].fail_user++;
                }
            }
        }
        Arrays.sort(stage);
        int[] answer = new int[stage.length];
        for(int i = 0; i < stage.length; i++) answer[i] = stage[i].index;
        return answer;
    }
}