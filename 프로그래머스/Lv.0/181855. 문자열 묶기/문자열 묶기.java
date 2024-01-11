class Solution {
    public int solution(String[] strArr) {
        int[] count = new int[31];
        for(int i = 0; i < strArr.length; i++) count[strArr[i].length()]++;
        int answer = 0;
        for(int i = 0; i < count.length; i++) answer = Math.max(answer, count[i]);
        return answer;
    }
}