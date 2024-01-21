class Solution {
    public int solution(int[] array) {
        int answer = 0;
        for(int i = 0; i < array.length; i++){
            int curr = array[i];
            while(curr > 0){
                if(curr % 10 == 7) answer++;
                curr /= 10;
            }
        }
        return answer;
    }
}