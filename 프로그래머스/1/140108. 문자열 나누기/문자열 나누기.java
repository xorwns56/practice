class Solution {
    public int solution(String s) {
        int answer = 0;
        char[] chars = s.toCharArray();
        char x = '\0';
        int[] count = new int[2];
        for(int i = 0; i < chars.length; i++){
            if(x == '\0') x = chars[i];
            if(chars[i] == x) count[0]++;
            else count[1]++;
            if(count[0] == count[1]){
                answer++;
                x = '\0';
                count[0] = count[1] = 0;
            }
        }
        return answer + (count[0] + count[1] > 0 ? 1 : 0);
    }
}