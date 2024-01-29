class Solution {
    public int solution(String[] babbling) {
        int answer = 0;
        for(int i = 0; i < babbling.length; i++){
            String s = babbling[i].replaceFirst("aya", " ").replaceFirst("ye", " ").replaceFirst("woo", " ").replaceFirst("ma", " ");
            if(s.trim().isEmpty()) answer++;
        }
        return answer;
    }
}