/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.awt.Desktop;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/** Explicitly opens only a local packet. Does not automate a website or inspect browser sessions. */
final class StudyBrowserLauncher {
    private StudyBrowserLauncher() { }

    static String open(Path handoff, String browser) throws IOException {
        Path dir = handoff.toRealPath(), page = dir.resolve("index.html");
        if (!Files.isRegularFile(dir.resolve("packet.yaml")) || !Files.isRegularFile(page) || Files.isSymbolicLink(page))
            throw new IOException("Expected a local handoff packet directory");
        StudyBrowserHandoff.verify(dir);
        if (!Set.of("auto", "chrome", "chromium", "firefox", "edge", "links", "none").contains(browser)) throw new IllegalArgumentException("Unknown browser selection");
        String fallback = "Open this local packet on an authorized workstation: " + page;
        if (browser.equals("none")) return fallback;
        if (browser.equals("auto") && Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE)) {
            try { Desktop.getDesktop().browse(page.toUri()); return "Requested the system-default browser for " + page + ". Browser rendering and submission are not confirmed."; }
            catch (IOException | UnsupportedOperationException | SecurityException unavailable) { /* Try an installed browser, never a remote service. */ }
        }
        String os = System.getProperty("os.name").toLowerCase(Locale.ROOT);
        boolean graphical = os.contains("windows") || os.contains("mac") || System.getenv("DISPLAY") != null || System.getenv("WAYLAND_DISPLAY") != null;
        List<String> choices = browser.equals("auto") ? graphical ? List.of("chrome", "chromium", "firefox", "edge", "links") : List.of("links") : List.of(browser);
        for (String choice : choices) for (Path candidate : candidates(choice, os, System.getenv())) {
            if (!Files.isRegularFile(candidate) || (!os.contains("windows") && !Files.isExecutable(candidate))) continue;
            if (choice.equals("links") && System.console() == null) return "Links requires an interactive terminal. " + fallback;
            try {
                // No shell, user-supplied command, remote URL, prompt argument or profile/cookie access.
                ProcessBuilder process = new ProcessBuilder(candidate.toAbsolutePath().toString(), page.toUri().toString());
                if (choice.equals("links")) process.inheritIO();
                else process.redirectOutput(ProcessBuilder.Redirect.DISCARD).redirectError(ProcessBuilder.Redirect.DISCARD);
                Process started = process.start();
                if (choice.equals("links")) {
                    try { started.waitFor(); }
                    catch (InterruptedException interrupted) { Thread.currentThread().interrupt(); return "Links viewer wait interrupted; no prompt submitted. " + fallback; }
                }
                return "Requested " + choice + " (" + candidate + ") for " + page + ". No prompt submitted."
                    + (choice.equals("links") ? " Links is a local packet viewer, not a supported client for modern AI sites." : "");
            } catch (IOException unavailable) { /* Continue bounded installed candidates. */ }
        }
        return "No usable " + browser + " browser found. " + fallback;
    }

    static List<Path> candidates(String browser, String os, Map<String, String> env) {
        List<String> names = switch (browser) {
            case "chrome" -> List.of("google-chrome", "google-chrome-stable", "chrome");
            case "chromium" -> List.of("chromium", "chromium-browser");
            case "firefox" -> List.of("firefox");
            case "edge" -> List.of("microsoft-edge", "msedge");
            case "links" -> List.of("links", "links2");
            default -> throw new IllegalArgumentException("Unknown browser");
        };
        List<Path> result = new ArrayList<>();
        boolean windows = os.contains("windows");
        String path = env.entrySet().stream().filter(e -> e.getKey().equalsIgnoreCase("PATH")).map(Map.Entry::getValue).findFirst().orElse("");
        for (String part : path.split(windows ? ";" : ":")) {
            if (part.isBlank()) continue;
            Path directory = Path.of(part.replaceAll("^\"|\"$", ""));
            if (!directory.isAbsolute()) continue;
            for (String name : names) result.add(directory.resolve(name + (windows ? ".exe" : "")));
        }
        if (windows) {
            String suffix = switch (browser) {
                case "chrome" -> "Google/Chrome/Application/chrome.exe";
                case "chromium" -> "Chromium/Application/chrome.exe";
                case "firefox" -> "Mozilla Firefox/firefox.exe";
                case "edge" -> "Microsoft/Edge/Application/msedge.exe";
                default -> null;
            };
            if (suffix != null) for (String key : List.of("ProgramFiles", "ProgramFiles(x86)", "LOCALAPPDATA"))
                if (env.containsKey(key)) result.add(Path.of(env.get(key)).resolve(suffix));
        } else if (os.contains("mac")) {
            String binary = switch (browser) {
                case "chrome" -> "Google Chrome.app/Contents/MacOS/Google Chrome";
                case "chromium" -> "Chromium.app/Contents/MacOS/Chromium";
                case "firefox" -> "Firefox.app/Contents/MacOS/firefox";
                case "edge" -> "Microsoft Edge.app/Contents/MacOS/Microsoft Edge";
                default -> null;
            };
            if (binary != null) result.add(Path.of("/Applications").resolve(binary));
        }
        return result;
    }
}
