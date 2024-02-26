import java.util.*;
class Solution {
    public int solution(int cacheSize, String[] cities) {
        int answer = 0;
        LinkedList<String> cache = new LinkedList<>();
        for(int i = 0; i < cities.length; i++){
            String city = cities[i].toLowerCase();
            int cache_index = cache.indexOf(city);
            answer += cache_index < 0 ? 5 : 1;
            if(cache_index >= 0) cache.remove(cache_index);
            cache.addFirst(city);
            if(cache.size() > cacheSize) cache.removeLast();
        }
        return answer;
    }
}