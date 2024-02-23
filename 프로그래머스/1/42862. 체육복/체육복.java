class Solution {
    public int solution(int n, int[] lost, int[] reserve) {
        int[] count = new int[n + 1];
        for(int i = 0; i < reserve.length; i++) count[reserve[i]] = 1;
        for(int i = 0; i < lost.length; i++) count[lost[i]]--;
        count[0] = -1;
        int answer = 0;
        for(int i = 1; i <= n; i++){
            if(count[i] >= 0) answer++;
            else{
                if(count[i - 1] > 0){
                    answer++;
                    count[i - 1]--;
                }else if(i + 1 <= n && count[i + 1] > 0){
                    answer++;
                    count[i + 1]--;
                }
            }
        }
        return answer;
    }
}