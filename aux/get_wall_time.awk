{
    if ($0 ~/wall clock/) {
        printf substr($0, index($0, "Elapsed (wall clock) time (h:mm:ss or m:ss): ") + length("Elapsed (wall clock) time (h:mm:ss or m:ss): "))
    }
}
