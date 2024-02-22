class Solution {
    public int solution(String dartResult) {
        char[] chars = dartResult.toCharArray();
        int answer = 0;
        int prev = 0;
        int i = 0;
        while(i < chars.length){
            int score = 0;
            while('0' <= chars[i] && chars[i] <= '9'){
                score = 10 * score + chars[i] - '0';
                i++;
            }
            int total_score = 0;
            while(chars[i]=='S' || chars[i]=='D' || chars[i]=='T' || chars[i]=='*' || chars[i]=='#'){
                switch(chars[i]){
                    case 'S': total_score = score; break;
                    case 'D': total_score = score * score; break;
                    case 'T': total_score = score * score * score; break;
                    case '*': total_score *= 2; prev *= 2; break;
                    case '#': total_score *= -1; break; 
                }
                i++;
                if(i >= chars.length) break;
            }
            answer += prev;
            prev = total_score;
        }
        return answer + prev;
    }
}