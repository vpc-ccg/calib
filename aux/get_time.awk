{
    if ($0 ~/User time/) {
        printf substr($0, index($0, "User time (seconds): ") + length("User time (seconds): "))
    }
}
