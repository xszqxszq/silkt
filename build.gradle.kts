plugins {
    kotlin("jvm") version "2.0.21"
    `java-library`
    `maven-publish`
    signing
    id("com.vanniktech.maven.publish") version "0.34.0"
}

java {
    withSourcesJar()
}

group = "xyz.xszq"
version = "1.0.0"

repositories {
    mavenCentral()
}

dependencies {
    testImplementation(kotlin("test"))
}

kotlin {
    jvmToolchain(22)
}

tasks.test {
    useJUnitPlatform()
}

mavenPublishing {
    publishToMavenCentral()
    signAllPublications()

    coordinates("xyz.xszq", "silkt", version as String)

    pom {
        name.set("silkt")
        description.set("Pure Kotlin SILK v3 encoder.")
        inceptionYear.set("2026")
        url.set("https://github.com/xszqxszq/silkt")
        licenses {
            license {
                name.set("The MIT License")
                url.set("https://opensource.org/license/mit")
            }
        }
        developers {
            developer {
                id.set("xszqxszq")
                name.set("xszqxszq")
                url.set("https://github.com/xszqxszq/")
            }
        }
        scm {
            url.set("https://github.com/xszqxszq/silkt/")
            connection.set("scm:git:git://github.com/xszqxszq/silkt.git")
            developerConnection.set("scm:git:ssh://git@github.com/xszqxszq/silkt.git")
        }
    }
}
