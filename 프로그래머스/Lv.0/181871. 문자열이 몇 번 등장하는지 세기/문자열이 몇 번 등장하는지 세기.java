class Solution {
    public int solution(String myString, String pat) {
        int count = 0;
        for(int i = 0; i < myString.length(); i++){
            if(pat.length() - 1 <= i){
                for(int j = 0; j < pat.length(); j++){
                    if(myString.charAt(i - (pat.length() - 1) + j) == pat.charAt(j)){
                        if(j == pat.length() - 1) count++;
                    }else break;
                }
            }
        }
        return count;
    }
}