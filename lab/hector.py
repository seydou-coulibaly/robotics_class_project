def check_nonsingularity(r5, commands) :
    q = 0
    length = len(commands)
    while q < length :
        del commands[q]
        if r5.actuate(commands[q]):
            length -= 1
        q += 1
    return commands
