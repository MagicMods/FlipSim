## TickLog Class Introduction (Sim/src/util/ticklog.js)

**Purpose:** To provide a simple time-based trigger mechanism for throttling high-frequency logging.

**Functionality:**

- An instance is created with a specific interval (`intervalMs`).
- The `update()` method, when called periodically, sets an internal boolean flag (`tick`) to `true` if the specified interval has passed since the last time it was set _and_ the flag is currently `false`.
- The `GetTick()` method returns the current state of the `tick` flag.
- The `ResetTick()` method sets the `tick` flag back to `false`.

**Usage Pattern:**

1.  Instantiate `TickLog` for each distinct log stream requiring throttling: `const myThrottler = new TickLog(intervalMs);`
2.  Call `myThrottler.update()` regularly within a loop or timer.
3.  In the code section performing the logging, check the flag: `if (myThrottler.GetTick()) { ... }`
4.  If the flag is `true`, perform the log operation.
5.  Immediately after logging, reset the flag: `myThrottler.ResetTick();`

This design places the responsibility of resetting the flag on the consuming code after the log action is taken. One instance is needed per log stream.
