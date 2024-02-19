class Solution {
    public int solution(String t, String p) {
        int answer = 0;
        char[] t_chars = t.toCharArray();
        char[] p_chars = p.toCharArray();
        for(int i = 0; i <= t_chars.length - p_chars.length; i++){
            for(int j = i; j < i + p_chars.length; j++){
                if(t_chars[j] > p_chars[j - i]) break;
                else if(t_chars[j] < p_chars[j - i] || j + 1 == i + p_chars.length){
                    answer++;
                    break;
                }
            }
        }
        return answer;
    }
}